// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoLazyNestedPassiveScalarFieldList.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sat January 9 2021
/// @brief    [\ref Enzo] Implementation of EnzoNestedPassiveScalarFieldList
///           class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

// define a shorthand alias type just for use in this f
typedef std::vector<std::string> str_vec_t;

//----------------------------------------------------------------------

std::shared_ptr<const std::vector<std::vector<std::string>>>
EnzoNestedPassiveScalarFieldList::build_nested_list_() noexcept
{

  std::shared_ptr<std::vector<str_vec_t>> nested_list =
    std::make_shared<std::vector<str_vec_t>>();

  str_vec_t groups = EnzoCenteredFieldRegistry::passive_scalar_group_names();

  // THIS FUNCTION NEEDS TO BE UPDATED TO PROPERLY HANDLE GROUPS OF FIELDS
  // THAT MUST HAVE MASS FRACTIONS THAT SUM TO ONE. The following assertion
  // error has been put into place to remind the developer of this.
  ASSERT("EnzoNestedPassiveScalarFieldList::build_nested_list_",
         "This function needs to be updated because additional groups of "
         "passive scalars were introduced.", groups.size() == 1)

  Grouping *grouping = cello::field_descr()->groups();

  nested_list->push_back(str_vec_t());

  for (std::string group_name : groups){
    int num_fields = grouping->size(group_name);
    for (int field_ind = 0; field_ind < num_fields; field_ind++){
      std::string field_name = grouping->item(group_name, field_ind);
      (*nested_list)[0].push_back(field_name);
    }
  }

  return std::const_pointer_cast<const std::vector<str_vec_t>>(nested_list);
}

//----------------------------------------------------------------------

void EnzoNestedPassiveScalarFieldList::pup(PUP::er &p) {
  p | initialized_;

  if (p.isUnpacking()){
    std::vector<str_vec_t>* temp_vals = new std::vector<str_vec_t>();
    p | *temp_vals;
    nested_names_ = std::const_pointer_cast<const std::vector<str_vec_t>>
      (std::shared_ptr<std::vector<str_vec_t>>(temp_vals));
  } else {
    std::vector<str_vec_t> temp_copy = *nested_names_;
    p | temp_copy;
  }
}
