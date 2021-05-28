// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoLazyNestedPassiveScalarFieldList.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sat January 9 2021
/// @brief    [\ref Enzo] Implementation of EnzoLazyPassiveScalarFieldList
///           class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

//----------------------------------------------------------------------

std::shared_ptr<const str_vec_t>
EnzoLazyPassiveScalarFieldList::build_passive_list_() noexcept
{
  std::shared_ptr<str_vec_t> passive_list = std::make_shared<str_vec_t>();
  str_vec_t groups = EnzoCenteredFieldRegistry::passive_scalar_group_names();
  Grouping *grouping = cello::field_descr()->groups();

  for (std::string group_name : groups){
    int num_fields = grouping->size(group_name);
    for (int field_ind = 0; field_ind < num_fields; field_ind++){
      std::string field_name = grouping->item(group_name, field_ind);
      passive_list->push_back(field_name);
    }
  }

  return std::const_pointer_cast<const str_vec_t>(passive_list);
}

//----------------------------------------------------------------------

void EnzoLazyPassiveScalarFieldList::pup(PUP::er &p) {
  p | initialized_;

  if (initialized_){
    if (p.isUnpacking()){
      str_vec_t* temp_vals = new std::vector<std::string>();
      p | *temp_vals;
      passive_names_ = std::const_pointer_cast<const str_vec_t>
	(std::shared_ptr<str_vec_t>(temp_vals));
    } else {
      str_vec_t temp_copy = *passive_names_;
      p | temp_copy;
    }
  }
}
