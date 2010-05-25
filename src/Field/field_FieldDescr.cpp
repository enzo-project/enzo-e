// $Id: field_FieldDescr.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the FieldDescr class

#include "error.hpp"
#include "field.hpp"

//----------------------------------------------------------------------

FieldDescr::FieldDescr () throw ()
  : alignment_(0),
    padding_(0),
    courant_(1),
    field_name_(),
    field_id_(),
    group_name_(),
    group_id_(),
    field_in_group_(),
    precision_(),
    centering_(),
    ghosts_(),
    min_value_(),
    max_value_(),
    min_action_(),
    max_action_()
{
}

//----------------------------------------------------------------------

int FieldDescr::field_count() const throw()
{
  return field_name_.size();
}

//----------------------------------------------------------------------

std::string FieldDescr::field_name(size_t id_field) const throw()
{ 
  return field_name_.at(id_field);
}

//----------------------------------------------------------------------

int FieldDescr::field_id(const std::string name) const throw()
{
  return field_id_.at(name);
}

//----------------------------------------------------------------------

int FieldDescr::group_count() const throw()
{
  return group_name_.size(); 
}

//----------------------------------------------------------------------

std::string FieldDescr::group_name(int id_group) const throw()
{
  return group_name_.at(id_group);
}

//----------------------------------------------------------------------

int FieldDescr::group_id(const std::string name) const throw()
{
  return group_id_.at(name);
}

//----------------------------------------------------------------------

bool FieldDescr::field_in_group(int id_field, int id_group) const throw()
{
  set_int_type t = field_in_group_[id_field];
  return t.find(id_group) != t.end();
}

//----------------------------------------------------------------------

int FieldDescr::alignment() const throw()
{
  return alignment_;
}

//----------------------------------------------------------------------

int FieldDescr::padding() const throw()
{
  return padding_;
}

//----------------------------------------------------------------------

void FieldDescr::centering
(
 int id_field,
 bool * cx, 
 bool * cy, 
 bool * cz
) const throw()
{
  *cx = centering_.at(id_field)[0];
  *cy = centering_.at(id_field)[1];
  *cz = centering_.at(id_field)[2];
}

//----------------------------------------------------------------------

void FieldDescr::ghosts
(
 int id_field,
 int * gx, 
 int * gy, 
 int * gz
) const throw()
{
  *gx = ghosts_.at(id_field)[0];
  *gy = ghosts_.at(id_field)[1];
  *gz = ghosts_.at(id_field)[2];
}

//----------------------------------------------------------------------

precision_type FieldDescr::precision(int id_field) const throw()
{
  return precision_[id_field];
}
