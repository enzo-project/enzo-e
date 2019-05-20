// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCenteredFieldRegistry.cpp
/// @author   Matthew W. Abruzzo (mwa2113@columbia.edu)
/// @date     Thurs March 28 2014
/// @brief    Implements the EnzoCenteredFieldRegistry class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

Grouping* EnzoCenteredFieldRegistry::build_grouping_
(const EnzoFieldConditions cond, const std::string leading_prefix,
 const std::string overlap_cat_prefix, FieldCat target_cat)
{
  Grouping *out = new Grouping;
  #define ENTRY(name, req_name, category, quantity_type)                      \
    if (cond.req_name && category_overlap_(category, target_cat)) {           \
      add_group_fields_(out, std::string( #name ),                            \
                        std::string( #quantity_type ), category,              \
                        leading_prefix, overlap_cat_prefix);		      \
    }
  FIELD_TABLE
  #undef ENTRY
  return out;
}

//----------------------------------------------------------------------

void EnzoCenteredFieldRegistry::add_group_fields_
(Grouping *grouping, const std::string group_name,
 const std::string quantity_type, FieldCat category,
 const std::string leading_prefix, const std::string overlap_cat_prefix)
{
  std::string prefix;
  if (category == FieldCat::overlap){
    prefix = leading_prefix + overlap_cat_prefix;
  } else {
    prefix = leading_prefix;
  }

  if (quantity_type == "SCALAR"){
    grouping->add(prefix + group_name, group_name);
  } else {
    grouping->add(prefix + group_name + std::string("_x"), group_name);
    grouping->add(prefix + group_name + std::string("_y"), group_name);
    grouping->add(prefix + group_name + std::string("_z"), group_name);
  }
}
//----------------------------------------------------------------------




std::vector<std::string> EnzoCenteredFieldRegistry::group_names_
 (const EnzoFieldConditions cond, bool include_passive, FieldCat target_cat)
{
  std::vector<std::string> out;
  #define ENTRY(name, req_name, category, quantity_type)                      \
    if (cond.req_name && category_overlap_(category, target_cat)) {           \
      out.push_back( std::string( #name ) );                                  \
    }
  FIELD_TABLE
  #undef ENTRY

  if (include_passive){
    std::vector<std::string> passive = passive_scalar_group_names();
    out.reserve(out.size() + passive.size());
    out.insert(out.end(), passive.begin(), passive.end());
  }
  return out;
}

//----------------------------------------------------------------------

#define PREPARE_LUT_SCALAR(name, req_name, in_category)                       \
  out.name = (cond.req_name && in_category) ? i++ : -1
#define PREPARE_LUT_VECTOR(name, req_name, in_category)                       \
  out.VECTOR_COMPONENT(name,_i) = (cond.req_name && in_category) ? i++ : -1;  \
  out.VECTOR_COMPONENT(name,_j) = (cond.req_name && in_category) ? i++ : -1;  \
  out.VECTOR_COMPONENT(name,_k) = (cond.req_name && in_category) ? i++ : -1

field_lut EnzoCenteredFieldRegistry::prepare_lut_
( const EnzoFieldConditions cond, int *nfields, FieldCat target_cat)
{
  field_lut out;
  int i = 0;
  #define ENTRY(name, req_name, category, quantity_type)                      \
    PREPARE_LUT_##quantity_type(name, req_name,                               \
                                category_overlap_(category, target_cat));
  FIELD_TABLE
  #undef ENTRY
  *nfields = i;
  return out;
}

//----------------------------------------------------------------------

#define USE_QUANTITY_SCALAR(lut, name) (lut.name > -1)
#define USE_QUANTITY_VECTOR(lut, name) (lut.VECTOR_COMPONENT(name,_i) > -1)

EFlt3DArray* EnzoCenteredFieldRegistry::load_array_of_fields
(Block *block, const field_lut lut, const int nfields, Grouping &grouping,
 const int dim)
{
  EFlt3DArray* arr = new EFlt3DArray[nfields];
  int cur_index = 0;
  EnzoPermutedCoordinates coord(dim);
  EnzoFieldArrayFactory array_factory(block);

  #define ENTRY(name, req_name, category, quantity_type)                      \
    if (USE_QUANTITY_##quantity_type(lut, name)) {		              \
      std::string group_name = std::string( #name );                          \
      std::string quantity_type = std::string( #quantity_type );              \
      cur_index = load_array_of_fields_(arr, cur_index, grouping, group_name, \
                                        quantity_type, dim, coord,            \
				        array_factory);		              \
  }
  FIELD_TABLE
  #undef ENTRY
  return arr;
}

//----------------------------------------------------------------------

int EnzoCenteredFieldRegistry::load_array_of_fields_
( EFlt3DArray* arr, int cur_index, Grouping &grouping, std::string group_name,
  std::string quantity_type, int dim, EnzoPermutedCoordinates coord,
  EnzoFieldArrayFactory array_factory)
{
  // Sanity Check:
  int group_size = grouping.size(group_name);

  ASSERT("load_fluid_fields_",
	 "Groups of fields don't have the correct number of fields.",
	 (quantity_type == "VECTOR" && group_size == 3) ||
	 (quantity_type == "SCALAR" && group_size == 1));

  if (quantity_type == "VECTOR"){
    arr[cur_index] = array_factory.reconstructed_field(grouping, group_name,
						       coord.i_axis(), dim);
    arr[cur_index+1] = array_factory.reconstructed_field(grouping,
							 group_name,
							 coord.j_axis(),
							 dim);
    arr[cur_index+2] = array_factory.reconstructed_field(grouping,
							 group_name,
							 coord.k_axis(),
							 dim);
    return cur_index+3;
  } else {

    arr[cur_index] = array_factory.reconstructed_field(grouping, group_name,
						       0, dim);
    return cur_index+1;
  }
}


