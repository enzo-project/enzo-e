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

bool EnzoCenteredFieldRegistry::is_actively_advected_vector_component
(std::string name, bool ijk_suffix) const noexcept
{
  const std::array<std::string,3> suffixes =  {ijk_suffix ? "_i" : "_x",
					       ijk_suffix ? "_j" : "_y",
					       ijk_suffix ? "_k" : "_z"};

  std::size_t length = name.length();
  if (length >= 2){
    std::string suffix = name.substr(length-2,2);
    if (std::find(suffixes.begin(), suffixes.end(), suffix) != suffixes.end()){
      std::string prefix = name.substr(0,length-2);
      bool is_vector, actively_advected;
      FieldCat category;
      bool success = quantity_properties(prefix, &is_vector, &category,
					 &actively_advected);
      return success && is_vector && actively_advected;
    }
  }
  return false;
}

//----------------------------------------------------------------------

std::string EnzoCenteredFieldRegistry::get_actively_advected_quantity_name
(std::string name, bool ijk_suffix) const noexcept
{
  std::string out = "";
  if (quantity_properties(name)){
    out = name;
  } else if (is_actively_advected_vector_component(name, true)){
    // current element is a VECTOR QUANTITY
    out = name.substr(0,name.length()-2);
  }
  return out;
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
