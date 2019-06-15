// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCenteredFieldRegistry.cpp
/// @author   Matthew W. Abruzzo (mwa2113@columbia.edu)
/// @date     Thurs March 28 2014
/// @brief    Implements the EnzoCenteredFieldRegistry class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"
#include <algorithm>
#include <tuple>
#include <map>


//----------------------------------------------------------------------

// Helper function used to check if a vector of strings contains an item
bool contains_item_(std::vector<std::string> names, std::string item)
{
  return (std::find(names.begin(), names.end(), item) != names.end());
}

//----------------------------------------------------------------------

typedef std::tuple<std::string, std::string, FieldCat, bool> FT_row;
// Runtime function to construct a dynamic version of FIELD_TABLE
// It may make sense down the road to cache this
std::vector<FT_row> get_dynamic_table_()
{
  std::vector<FT_row> out;
  #define ENTRY(name, math_type, category, if_advection)                      \
    if (std::string( #if_advection) == "T"){                                  \
      out.push_back(std::make_tuple( #row , #math_type, category, true));     \
    } else {                                                                  \
      out.push_back(std::make_tuple( #row , #math_type, category, false));    \
    }
  FIELD_TABLE
  #undef ENTRY
}

//----------------------------------------------------------------------

std::vector<std::string> EnzoCenteredFieldRegistry::get_registered_quantities()
{
  std::vector<FT_row> table = get_dynamic_table_();
  std::vector<std::string> out;
  for(auto const& row: table) {
    out.push_back(std::get<0>(row));
  }
  return out;
}

//----------------------------------------------------------------------

std::vector<std::string> EnzoCenteredFieldRegistry::get_registered_fields()
{
  std::vector<FT_row> table = get_dynamic_table_();
  std::vector<std::string> out;
  for(auto const& row: table) {
    std::string quantity_name = std::get<0>(row);
    if (std::get<1>(row) == "SCALAR"){
      out.push_back(quantity_name);
    } else{
      out.push_back(quantity_name + std::string("_x"));
      out.push_back(quantity_name + std::string("_y"));
      out.push_back(quantity_name + std::string("_z"));
    }
  }
  return out;
}

//----------------------------------------------------------------------

void EnzoCenteredFieldRegistry::check_known_quantity_names
(const std::vector<std::string> names)
{
  std::vector<std::string> known_names = get_registered_quantities();
  for(auto const& name: names) {
    if (!contains_item_(known_names, name)){
      ERROR1("EnzoCenteredFieldRegistry::check_known_names",
	     "%s is not a registered name", name.c_str());
    }
  }
}

//----------------------------------------------------------------------

Grouping* EnzoCenteredFieldRegistry::build_grouping
(const std::vector<std::string> quantity_names,
 const std::string leading_prefix = "")
{
  check_known_quantity_names(quantity_names);
  Grouping *out = new Grouping;
  #define ENTRY(name, math_type, category, if_advection)                      \
    if (contains_item_(quantity_names, std::string(name))){                   \
      add_group_fields_(out, std::string( #name ), std::string( #math_type ), \
			leading_prefix);                                      \
    }
  FIELD_TABLE
  #undef ENTRY
  return out;
}

//----------------------------------------------------------------------

void EnzoCenteredFieldRegistry::add_group_fields_
(Grouping *grouping, const std::string group_name,
 const std::string quantity_type, FieldCat category,
 const std::string leading_prefix)
{

  if (quantity_type == "SCALAR"){
    grouping->add(leading_prefix + group_name, group_name);
  } else {
    grouping->add(leading_prefix + group_name + std::string("_x"), group_name);
    grouping->add(leading_prefix + group_name + std::string("_y"), group_name);
    grouping->add(leading_prefix + group_name + std::string("_z"), group_name);
  }
}

//----------------------------------------------------------------------

void EnzoCenteredFieldRegistry::prepare_lut_
   const std::vector<std::string> quantity_names,
   int &conserved_start, int &conserved_stop, int &specific_start,
   int &specific_stop, int &other_start, int &other_stop, int &nfields
   const std::vector<std::string> flagged_quantities)
{
  // check that all entries flagged_quantities are also in quantity_names
  for(auto const& name: flagged_quantities) {
    if (!contains_item_(quantity_names, name)){
      ERROR1("EnzoCenteredFieldRegistry::categorize_quantities_",
	     "%s is listed in flagged_quantites but not in quantity_names",
	     name.c_str());
    }
  }

  // Main goal here is determine number conserved/specific/other fields
  // that are flagged/unflagged (track unflagged at 0 and flagged at 1)
  int num_conserved_fields[2] = {0, 0};
  int num_specific_fields[2] = {0, 0};
  int num_other_fields[2] = {0, 0};

  // We are going to do this without macros for simplicity
  // This is not necessarily the most efficient way to do things
  std::vector<FT_row> table = get_dynamic_table_();

  for(auto const& row: table) {
    std::string name, math_type;
    FieldCat category;
    bool if_advec;
    std::tie(name, math_type, category, if_advec) = row;

    if (!contains_item_(quantity_names, name)){
      continue;
    } else if (!if_advec) {
      ERROR1("EnzoCenteredFieldRegistry::categorize_quantities_",
	     ("quantity_names contains the field %s which is not an advection "
	      "related quantity."), name.c_str());
    }

    // If the quantity is flagged then update the value at index 1
    int update_ind = (contains_item_(flagged_quantities,name)) ? 1 : 0;

    // How much to increment the number of fields by
    int increment = (math_type == "VECTOR") 3 ? 1;
    if (category == FieldCat::conserved){
      num_conserved_fields[update_ind] += increment;
    } else if (category == FieldCat::specific){
      num_specific_fields[update_ind] += increment;
    } else {
      num_other_fields[update_ind] += increment;
    }
  }

  // finally let's set the returned information
  conserved_start = 0;
  conserved_stop  = conserved_start + num_conserved_fields[0];
  specific_start  = conserved_stop  + num_conserved_fields[1];
  specific_stop   = specific_start  + num_specific_fields[0];
  other_start     = specific_stop   + num_specific_fields[1];
  other_stop      = other_start     + num_other_fields[0];

  nfields         = other_stop      + num_other_fields[1];
}


//----------------------------------------------------------------------
/// @def      PREPARE_LUT
/// @brief    Macros that help initialize the values of the
///           EnzoAdvectionFieldLUT class
#define PREPARE_LUT_SCALAR(out, name, quantity_names, i)                      \
  out.name = (contains_item_(quantity_names,#name)) ? (*i)++ : -1
#define PREPARE_LUT_VECTOR(out, name, quantity_names, i)                      \
  out.COMBINE(name,_i) = (contains_item_(quantity_names,#name)) ? (*i)++ : -1;\
  out.COMBINE(name,_j) = (contains_item_(quantity_names,#name)) ? (*i)++ : -1;\
  out.COMBINE(name,_k) = (contains_item_(quantity_names,#name)) ? (*i)++ : -1


//----------------------------------------------------------------------
/// @def      PREPARE_ADVECTION_LUT
/// @brief    Macro that helps initialize a member of a member of the
///           EnzoAdvectionFieldLUT class
/// @param out the struct to be modified
/// @param name the name of the quantity to add (this should be a macro
///  argument taken from FIELD_TABLE)
/// @param math_type identification of the quantity as SCALAR or VECTOR (this
///  should be a macro argument taken from FIELD_TABLE)
/// @param quantity_names a vector of strings of quantity names that are
///  to have valid indices (fields not included in the table are set to -1
/// @param i pointer to the value that will indicate the index of the member
///  and that needs to be incremented. If the stringized version of name is not
///  contained by quantity_names then nothing is done. This value will be
///  different based whether the quantity is conserved, specific, or other and
///  whether it's flagged
#define PREPARE_ADVECTION_LUT_T(out, name, math_type, quantity_names, i)      \
  PREPARE_LUT_##math_type(out, name, quantity_names, ind)
#define PREPARE_ADVECTION_LUT_F(out, name, math_type, quantity_names, i)      \
  /* ... */


//----------------------------------------------------------------------
EnzoAdvectionFieldLUT EnzoCenteredFieldRegistry::prepare_advection_lut
  (const std::vector<std::string> quantity_names,
   int &conserved_start, int &conserved_stop, int &specific_start,
   int &specific_stop, int &other_start, int &other_stop, int &nfields,
   const std::vector<std::string> flagged_quantities)
{

  // First determin the values of conserved_start, conserved_stop, ..., nfields
  prepare_lut_(quantity_names, conserved_start, conserved_stop, specific_start,
	       specific_stop, other_start, other_stop, nfields,
	       flagged_quantities);

  // define the following indices to help setup the values of the look up table
  int conserved_index      = conserved_start;
  int flag_conserved_index = conserved_stop;
  int specific_index       = specific_start;
  int flag_specific_index  = specific_stop;
  int other_index          = other_start;
  int flag_other_index     = other_stop;

  EnzoAdvectionFieldLUT out;

  int *cur_ind;
  bool is_flagged;

  #define ENTRY(name, math_type, category, if_advection)                      \
    /* first we need to set cur_ind to the appropriate value */               \
    is_flagged = contains_item_(flagged_quantities, #name);                   \
    if (category == FieldCat::conserved) {                                    \
      cur_ind = (is_flagged) ? &flag_conserved_index : &conserved_index;      \
    } else if (category == FieldCat::specific) {                              \
      cur_ind = (is_flagged) ? &flag_specific_index : &specific_index;        \
    } else { /* category = FieldCat::other */                                 \
      cur_ind = (is_flagged) ? &flag_other_index : &other_index;              \
    }                                                                         \
                                                                              \
    /* Now use cur_ind to set value of out's member(s) (if applicable) */     \
    PREPARE_ADVECTION_LUT_##if_advection(out, name, math_type,                \
					 quantity_names, cur_ind);
  FIELD_TABLE
  #undef ENTRY

  return out;
}

//----------------------------------------------------------------------

#define USE_QUANTITY_SCALAR(lut, name) (lut.name > -1)
#define USE_QUANTITY_VECTOR(lut, name) (lut.COMBINE(name,_i) > -1)

#define LOAD_FIELDS_SCALAR(name, lut, arr, grouping, dim, coord, recon, arr_f)\
  load_array_of_fields_(arr, lut.name, grouping, name, "SCALAR", dim, 0       \
                        recon, arr_f)

#define LOAD_FIELDS_VECTOR(name, lut, arr, grouping, dim, coord, recon, arr_f)\
  load_array_of_fields_(arr, lut.COMBINE(name,_i), grouping, name, "VECTOR",  \
                        dim, coord.i_axis(), recon, arr_f);
  load_array_of_fields_(arr, lut.COMBINE(name,_j), grouping, name, "VECTOR",  \
                        dim, coord.j_axis(), recon, arr_f);
  load_array_of_fields_(arr, lut.COMBINE(name,_k), grouping, name, "VECTOR",  \
                        dim, coord.k_axis(), recon, arr_f)

// mt = math type
#define LOAD_FIELDS_T(name, mt, lut, arr, grouping, dim, coord, recon, arr_f) \
  if (USE_QUANTITY_##mt){                                                     \
    LOAD_FIELDS_##mt(name, lut, arr, grouping, dim, coord, recon, arr_f)l;    \
  }

#define LOAD_FIELDS_F(name, mt, lut, arr, grouping, dim, coord, recon, arr_f) \
  /* ... */

// This function should probably be refactor to not use macros at all. Instead
// we could just use the unary_advec_struct_for_each_ function defined in the
// header file

EFlt3DArray* EnzoCenteredFieldRegistry::load_array_of_fields
(Block *block, const EnzoAdvectionFieldLUT lut, const int nfields,
 Grouping &grouping, const int dim, const int stale_depth)
{
  EFlt3DArray* arr = new EFlt3DArray[nfields];
  EnzoPermutedCoordinates coord(dim);
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  // Set this to false if the loaded field is not reconstructed
  bool reconstructed = (dim != -1);

  #define ENTRY(name, math_type, category, if_advection)                      \
    LOAD_FIELDS_T(name, math_type, lut, arr, grouping, dim, coord,            \
                  reconstructed, array_factory)
  FIELD_TABLE
  #undef ENTRY
  return arr;
}

//----------------------------------------------------------------------

void EnzoCenteredFieldRegistry::load_array_of_fields_
( EFlt3DArray* arr, int index, Grouping &grouping, std::string group_name,
  std::string quantity_type, int dim, int group_ind, bool reconstructed,
  EnzoFieldArrayFactory array_factory)
{
  // Sanity Check:
  int group_size = grouping.size(group_name);
  ASSERT("load_array_of_fields_",
	 "Groups of fields don't have the correct number of fields.",
	 (quantity_type == "VECTOR" && group_size == 3) ||
	 (quantity_type == "SCALAR" && group_size == 1));

  if (reconstructed){
    arr[index] = array_factory.reconstructed_field(grouping, group_name,
						   group_ind, dim);
  } else {
    arr[index] = array_factory.from_grouping(grouping, group_name, group_ind);
  }
}


