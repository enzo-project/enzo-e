// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTempFieldUtils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Tues April 28 2020
/// @brief    [\ref Enzo] Implementation of functions in the EnzoTempFieldUtils
///           namespace

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoTempFieldUtils::prep_reused_temp_field
(Field &field, std::string field_name, int cx, int cy, int cz)
{
  FieldDescr * field_descr = field.field_descr();
  int id;
  if (field.is_field(field_name)){
    // the temporary field was added in a previous time-step
    id = field.field_id(field_name);
  } else {
    // Reserve a temporary field
    id = field_descr->insert_temporary(field_name);
    field_descr->set_centering(id,cx,cy,cz);
  }

  // allocate temporary field
  // this last check allows us to make temporary field "permanent" for debugging
  // purposes
  if (!field.is_permanent(id)){
    field.allocate_temporary(id);
  }
}

//----------------------------------------------------------------------

void EnzoTempFieldUtils::deallocate_grouping_fields
(Field &field, const std::vector<std::string> &group_names, Grouping &grouping)
{
  for (unsigned int i=0;i<group_names.size();i++){
    std::string group_name = group_names[i];
    int num_fields = grouping.size(group_name);
    for (int j=0;j<num_fields;j++){
      // Determine field_name and id
      std::string field_name = (grouping.item(group_name,j));
      int id = field.field_id(field_name);
      // if it's a temporary field, deallocate it
      if (!field.is_permanent(id)) {
	field.deallocate_temporary(id);
      }
    }
  }
}
