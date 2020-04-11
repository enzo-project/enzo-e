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

void print_lut(const EnzoAdvectionFieldLUT lut)
{

  std::string out = std::string("{\n");

  // Define the lambda function to collect values and names
  auto func = [&out](std::string field_name, int val)
    {
      std::string str_val = std::to_string(val);
      out = (out + std::string("  ") + field_name + std::string(" = ") +
	     str_val + std::string(";\n"));
    };

  // apply lambda function
  unary_advec_struct_for_each_(lut, func);
  CkPrintf("%s}\n", out.c_str());
}

//----------------------------------------------------------------------

void print_looked_up_vals(const EnzoAdvectionFieldLUT lut,
			  const enzo_float* array)
{

  std::string out = std::string("{\n");

  // Define the lambda function to collect values and names
  auto func = [&out, array](std::string field_name, int index)
    {
      std::string str_val = std::string("N/A");
      if (index != -1){
	char temp[25];
	sprintf(temp, "%.15e", array[index]);
	str_val = std::string(temp);
      }
      out = (out + std::string("  ") + field_name + std::string(" = ") +
	     str_val + std::string(";\n"));
    };

  unary_advec_struct_for_each_(lut, func);
  CkPrintf("%s}\n", out.c_str());
}

//----------------------------------------------------------------------

// Helper function used to check if a vector of strings contains an item
bool contains_item_(std::vector<std::string> names, std::string item)
{
  return (std::find(names.begin(), names.end(), item) != names.end());
}

//----------------------------------------------------------------------

// Runtime function to construct a dynamic version of FIELD_TABLE
// It may make sense down the road to cache this
std::map<std::string,FT_row> get_dynamic_table_(std::vector<std::string> &keys)
{
  std::map<std::string,FT_row> out;
  #define ENTRY(name, math_type, category, if_advection)                      \
    if (std::string( #if_advection) == "T"){                                  \
      out[#name] = std::make_tuple( #math_type , category,  true);            \
      keys.push_back(#name);                                                  \
    } else {                                                                  \
      out[#name] = std::make_tuple( #math_type , category, false);            \
      keys.push_back(#name);                                                  \
    }
  FIELD_TABLE
  #undef ENTRY
  return out;
}

//----------------------------------------------------------------------

EnzoCenteredFieldRegistry::EnzoCenteredFieldRegistry()
{
  field_table_ = get_dynamic_table_(table_keys_);
}

//----------------------------------------------------------------------

std::vector<std::string> EnzoCenteredFieldRegistry::get_registered_fields()
  const
{
  std::vector<std::string> out;
  for(auto const& key: table_keys_) {
    std::string quantity_name = key;
    if (std::get<0>(field_table_.at(key)) == "SCALAR"){
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
(const std::vector<std::string> names) const
{
  for(auto const& name: names) {
    if (!contains_item_(table_keys_, name)){
      ERROR1("EnzoCenteredFieldRegistry::check_known_names",
	     "%s is not a registered name", name.c_str());
    }
  }
}

//----------------------------------------------------------------------

bool EnzoCenteredFieldRegistry::quantity_properties
(std::string name, bool *vector_quantity, FieldCat *category,
 bool *actively_advected) const
{
  if (contains_item_(table_keys_, name)){
    FT_row row = field_table_.at(name);
    if (vector_quantity){
      *vector_quantity = ("VECTOR" == std::get<0>(row));
    }
    if (category){
      *category = std::get<1>(row);
    }
    if (actively_advected){
      *actively_advected = std::get<2>(row);
    }
    return true;
  } else {
    return false;
  }
}

//----------------------------------------------------------------------

Grouping* EnzoCenteredFieldRegistry::build_grouping
(const std::vector<std::string> quantity_names,
 const std::string leading_prefix) const
{
  check_known_quantity_names(quantity_names);
  Grouping *out = new Grouping;
  #define ENTRY(name, math_type, category, if_advection)                      \
    if (contains_item_(quantity_names, std::string(#name))){                  \
      add_group_fields_(out, std::string( #name ), std::string( #math_type ), \
			category, leading_prefix);                            \
    }
  FIELD_TABLE
  #undef ENTRY
  return out;
}

//----------------------------------------------------------------------

void EnzoCenteredFieldRegistry::add_group_fields_
(Grouping *grouping, const std::string group_name,
 const std::string quantity_type, FieldCat category,
 const std::string leading_prefix) const
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
(const std::vector<std::string> quantity_names,
 int &conserved_start, int &conserved_stop, int &specific_start,
 int &specific_stop, int &other_start, int &other_stop, int &nfields,
 const std::vector<std::string> flagged_quantities) const
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
  std::vector<std::string> keys = table_keys_;
  std::map<std::string, FT_row> table = field_table_;

  for(auto const& name: keys) {
    std::string math_type;
    FieldCat category;
    bool if_advec;
    std::tie(math_type, category, if_advec) = table[name];

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
    int increment = (math_type == "VECTOR") ? 3 : 1;
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


// determines the quantity name corresponding to a member of
// EnzoAdvectionFieldLUT
void parse_member_name_(std::string member_name, EnzoPermutedCoordinates coord,
			std::string &group_name, int &group_index)
{
  group_name = member_name;
  group_index = 0;
  std::size_t length = member_name.length();
  if (length >= 2){
    std::string suffix = member_name.substr(length-2,2);
    if (suffix == std::string("_i")){
      group_index = coord.i_axis();
    } else if (suffix == std::string("_j")){
      group_index = coord.j_axis();
    } else if (suffix == std::string("_k")){
      group_index = coord.k_axis();
    } else {
      return;
    }

    group_name = member_name.substr(0,length-2);
  }
}

//----------------------------------------------------------------------

// helper function that determines the quantity name (or group name) that
// corresponds to the name of a member of EnzoAdvectionFieldLUT
std::string determine_quantity_name_(std::string member_name)
{
  EnzoPermutedCoordinates temp_coord(0);
  int temp_index;
  std::string quantity_name;
  parse_member_name_(member_name, temp_coord, quantity_name, temp_index);
  return quantity_name;
}


//----------------------------------------------------------------------

EnzoAdvectionFieldLUT EnzoCenteredFieldRegistry::prepare_advection_lut
  (const std::vector<std::string> quantity_names,
   int &conserved_start, int &conserved_stop, int &specific_start,
   int &specific_stop, int &other_start, int &other_stop, int &nfields,
   const std::vector<std::string> flagged_quantities) const
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

  std::vector<std::string> keys = table_keys_;
  std::map<std::string,FT_row> table = field_table_;

  // define a lambda function to use while iterating over the entries of out
  auto fn = [quantity_names, flagged_quantities, table,
	     &conserved_index, &flag_conserved_index,
	     &specific_index, &flag_specific_index,
	     &other_index, &flag_other_index](std::string name, int &val)
    {
      std::string quantity_name = determine_quantity_name_(name);

      if (!contains_item_(quantity_names, quantity_name)){
	// if the quantity name is not within quantity_names set value to -1
	val = -1;
      } else {
	bool is_flagged = contains_item_(flagged_quantities, quantity_name);
	FieldCat category = std::get<1>(table.at(quantity_name));
	if (category == FieldCat::conserved) {
	  val = (is_flagged) ? flag_conserved_index++ : conserved_index++;
	} else if (category == FieldCat::specific) {
	  val = (is_flagged) ?  flag_specific_index++ :  specific_index++;
	} else {
	  val = (is_flagged) ?     flag_other_index++ :     other_index++;
	}
      }
    };

  unary_advec_struct_for_each_(out, fn);
  return out;
}

//----------------------------------------------------------------------

EFlt3DArray* EnzoCenteredFieldRegistry::load_array_of_fields
(Block *block, const EnzoAdvectionFieldLUT lut, const int nfields,
 Grouping &grouping, const int dim, const int stale_depth) const
{
  EFlt3DArray* arr = new EFlt3DArray[nfields];
  // in the case where we don't have reconstructed values (dim = -1) we assume
  // that the that i-axis is aligned with the x-axis
  EnzoPermutedCoordinates coord( (dim == -1) ? 0 : dim);
  EnzoFieldArrayFactory array_factory(block, stale_depth);

  // define a lambda function to execute for every member of lut
  auto fn = [arr, coord, dim, &array_factory, &grouping](std::string name,
							int index)
    {
      // name is the name of a given member of the lut
      // index is the value associated with the member
      
      if (index != -1){
	int group_index;
	std::string group_name;
	parse_member_name_(name, coord, group_name, group_index);

	// Sanity Check:
	std::string quantity_type = (group_name == name) ? "SCALAR" : "VECTOR";
	int group_size = grouping.size(group_name);

	int expected_size = (quantity_type == "VECTOR") ? 3 : 1;
	ASSERT3("load_array_of_fields_",
	        "The \"%s\" group holds %d field(s). It should hold %d.",
		group_name.c_str(), group_size, expected_size,
		group_size == expected_size);

	if (dim != -1){
	  arr[index] = array_factory.reconstructed_field(grouping, group_name,
							 group_index, dim);
	} else {
	  arr[index] = array_factory.from_grouping(grouping, group_name,
						   group_index);
	}

      }
    };

  unary_advec_struct_for_each_(lut, fn);

  return arr;
}
