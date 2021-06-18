// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCenteredFieldRegistry.cpp
/// @author   Matthew W. Abruzzo (mwa2113@columbia.edu)
/// @date     Thurs March 28 2014
/// @brief    Implements the EnzoCenteredFieldRegistry class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"
#include <algorithm>

//----------------------------------------------------------------------
#define INIT_ROW_T_SCALAR(name, category) {#name , {false, category, true }},
#define INIT_ROW_T_VECTOR(name, category) {#name , {true,  category, true }},
#define INIT_ROW_F_SCALAR(name, category) {#name , {false, category, false}},
#define INIT_ROW_F_VECTOR(name, category) {#name , {true,  category, false}},

// initialize the static const map representing the field table
const ft_map EnzoCenteredFieldRegistry::field_table_ = {
    #define ENTRY(name, math_type, category, if_advection) \
    INIT_ROW_ ## if_advection ## _ ## math_type (name, category)
    FIELD_TABLE
    #undef ENTRY
};

//----------------------------------------------------------------------

std::vector<std::string> EnzoCenteredFieldRegistry::get_registered_quantities
(bool enumerate_components)
{
  std::vector<std::string> out;
  for (auto const &key_item_pair : EnzoCenteredFieldRegistry::field_table_) {
    if (enumerate_components && (key_item_pair.second).vector_quantity){
      out.push_back(key_item_pair.first + std::string("_x"));
      out.push_back(key_item_pair.first + std::string("_y"));
      out.push_back(key_item_pair.first + std::string("_z"));
    } else {
      out.push_back(key_item_pair.first);
    }
  }
  return out;
}

//----------------------------------------------------------------------

bool EnzoCenteredFieldRegistry::quantity_properties
(const std::string &name, bool *vector_quantity, FieldCat *category,
 bool *actively_advected) noexcept
{
  auto search = EnzoCenteredFieldRegistry::field_table_.find(name);
  if (search == EnzoCenteredFieldRegistry::field_table_.cend()){
    return false;
  } else {
    const ft_row &row = search->second;
    if (vector_quantity){ *vector_quantity = row.vector_quantity; }
    if (category){ *category = row.category; }
    if (actively_advected){ *actively_advected = row.actively_advected; }
    return true;
  }
}

//----------------------------------------------------------------------

bool is_actively_advected_vector_component_
(std::string name, bool ijk_suffix) noexcept
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
      bool success = EnzoCenteredFieldRegistry::quantity_properties
        (prefix, &is_vector, &category, &actively_advected);
      return success && is_vector && actively_advected;
    }
  }
  return false;
}

//----------------------------------------------------------------------

std::string EnzoCenteredFieldRegistry::get_actively_advected_quantity_name
(std::string name, bool ijk_suffix) noexcept
{
  std::string out = "";
  if (EnzoCenteredFieldRegistry::quantity_properties(name)){
    out = name;
  } else if (is_actively_advected_vector_component_(name, true)){
    // current element is a VECTOR QUANTITY
    out = name.substr(0,name.length()-2);
  }
  return out;
}
