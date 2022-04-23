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
(const bool enumerate_components, const bool exclude_non_active_advection)
{
  std::vector<std::string> out;
  for (auto const &key_item_pair : EnzoCenteredFieldRegistry::field_table_) {
    if (exclude_non_active_advection &&
        (!(key_item_pair.second).actively_advected)){
      continue;
    }
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

namespace{
  /// Identifies the vector component (if any) associated with a name.
  ///
  /// For all names consisting of 3 or more characters, this checks whether the
  /// last 2-characters in the string match a list of known suffixes (the
  /// list's contents depend on the `ijk_suffix` parameter). Upon success, this
  /// returns the character corresponding to the suffix. In all other cases,
  /// this returns the null character, `'\0'`.
  ///
  /// @param name The name that will be checked
  /// @param ijk_suffix When true, the suffixes are `{'_i', '_j', '_k'}`, and
  ///     the function returns `'i'`, `'j'`, or `'k'` upon a successful match.
  ///     Otherwise, the suffixes are `{'_x', '_y', '_z'}`, and the function
  ///     returns `'x'`, `'y'`, or `'z'` upon a successful match.
  char try_get_vector_component(const std::string& name,
                                bool ijk_suffix) noexcept
  {
    const std::array<char,3> ref = {ijk_suffix ? 'i' : 'x',
                                    ijk_suffix ? 'j' : 'y',
                                    ijk_suffix ? 'k' : 'z'};
    std::size_t len = name.length();
    if ((len >= 3) && (name[len - 2] == '_')){
        if (std::find(ref.begin(), ref.end(), name[len - 1]) != ref.end()){
            return name[len - 1];
        }
    }
    return '\0';
  }
}

//----------------------------------------------------------------------

std::string EnzoCenteredFieldRegistry::get_actively_advected_quantity_name
(const std::string& name, bool ijk_suffix, char* vec_component_rslt) noexcept
{

  // check if name matches the name of a scalar quantity)
  {
    bool is_vector, actively_advected;
    FieldCat category;
    bool success = EnzoCenteredFieldRegistry::quantity_properties
      (name, &is_vector, &category, &actively_advected);

    if (success && (!is_vector) && actively_advected) {
      if (vec_component_rslt != nullptr) { *vec_component_rslt = '\0'; }
      return name;
    }
  }

  // next, check if name corresponds to a component of a vector quantity
  const char vec_comp =  try_get_vector_component(name, ijk_suffix);
  if (vec_comp != '\0'){
    std::string candidate = name.substr(0, name.length()-2);

    bool is_vector, actively_advected;
    FieldCat category;
    bool success = EnzoCenteredFieldRegistry::quantity_properties
      (candidate, &is_vector, &category, &actively_advected);

    if (success && is_vector && actively_advected) {
      if (vec_component_rslt != nullptr) { *vec_component_rslt = vec_comp; }
      return candidate;
    }
  }

  // return empty string to indicate failure
  if (vec_component_rslt != nullptr) { *vec_component_rslt = '\0'; }
  return "";
}
