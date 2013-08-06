// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-14
/// @brief    Implementation of the FieldDescr class

#include "cello.hpp"

#include "field.hpp"

//----------------------------------------------------------------------

FieldDescr::FieldDescr () 
  throw ()
  : field_name_(),
    field_id_(),
    group_name_(),
    group_id_(),
    // field_in_group_(),
    alignment_(1),
    padding_(0),
    courant_(1),
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

FieldDescr::~FieldDescr() 
  throw()
{
  for (size_t i=0; i<centering_.size(); i++) {
    delete [] centering_[i];
    centering_[i] = 0;
  }
  for (size_t i=0; i<ghosts_.size();    i++) {
    delete [] ghosts_[i];
    ghosts_[i] = 0;
  }
}

//----------------------------------------------------------------------

FieldDescr::FieldDescr(const FieldDescr & field_descr) throw() 
{
  copy_(field_descr);
}

//----------------------------------------------------------------------

FieldDescr & FieldDescr::operator= (const FieldDescr & field_descr) 
  throw()
{
  copy_(field_descr);

  return *this;
}

//----------------------------------------------------------------------

int FieldDescr::field_count() const 
  throw()
{
  return field_name_.size();
}

//----------------------------------------------------------------------

std::string FieldDescr::field_name(size_t id_field) const 
  throw(std::out_of_range)
{ 
  return field_name_.at(id_field);
}

//----------------------------------------------------------------------

bool FieldDescr::is_field(const std::string & name) const 
  throw()
{ 
  return (field_id_.find(name) != field_id_.end());
}

//----------------------------------------------------------------------

int FieldDescr::field_id(const std::string & name) const
  throw()
{
  //  return field_id_[name]; // ERROR IN PGI ON GORDON 11.9-0 64-bit
  std::map<const std::string,int>::const_iterator it;
  it=field_id_.find(name);
  if (it != field_id_.end()) {
    return it->second;
  } else {
    return -1;
  }
}

//----------------------------------------------------------------------

int FieldDescr::group_count() const 
  throw()
{
  return group_name_.size(); 
}

//----------------------------------------------------------------------

std::string FieldDescr::group_name(int id_group) const 
  throw(std::out_of_range)
{
  return group_name_[id_group];
}

//----------------------------------------------------------------------

bool FieldDescr::is_group(const std::string & name) const
  throw()
{ 
  return (group_id_.find(name) != group_id_.end());
}

//----------------------------------------------------------------------

int FieldDescr::group_id(const std::string & name) const
  throw()
{
  std::map<const std::string,int>::const_iterator it;
  it=group_id_.find(name);
  if (it != group_id_.end()) {
    return it->second;
  } else {
    return -1;
  }
  //  return group_id_[name];
}

//----------------------------------------------------------------------

// bool FieldDescr::field_in_group(int id_field, int id_group) const 
//   throw(std::out_of_range)
// {
//   int_set_type t = field_in_group_.at(id_field);
//   return t.find(id_group) != t.end();
// }

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

double FieldDescr::courant() const throw()
{
  return courant_;
}

//----------------------------------------------------------------------

void FieldDescr::centering
(
 int id_field,
 bool * cx, 
 bool * cy, 
 bool * cz
 ) const throw(std::out_of_range)
{
  if (cx) (*cx) = centering_.at(id_field)[0];
  if (cy) (*cy) = centering_.at(id_field)[1];
  if (cz) (*cz) = centering_.at(id_field)[2];
}

//----------------------------------------------------------------------

void FieldDescr::ghosts
(
 int id_field,
 int * gx, 
 int * gy, 
 int * gz
 ) const throw(std::out_of_range)
{
  if (gx) (*gx) = ghosts_.at(id_field)[0];
  if (gy) (*gy) = ghosts_.at(id_field)[1];
  if (gz) (*gz) = ghosts_.at(id_field)[2];
}

//----------------------------------------------------------------------

precision_type FieldDescr::precision(int id_field) const 
  throw(std::out_of_range)
{
  return precision_.at(id_field);
}

//----------------------------------------------------------------------

double FieldDescr::minimum_value(int id_field) const 
  throw(std::out_of_range)
{
  return min_value_.at(id_field);
}

//----------------------------------------------------------------------

field_action_type FieldDescr::minimum_action(int id_field) const 
  throw(std::out_of_range)
{
  return min_action_.at(id_field);
}

//----------------------------------------------------------------------

double FieldDescr::maximum_value(int id_field) const 
  throw(std::out_of_range)
{
  return max_value_.at(id_field);
}

//----------------------------------------------------------------------

field_action_type FieldDescr::maximum_action(int id_field) const 
  throw(std::out_of_range)
{
  return max_action_.at(id_field);
}

//----------------------------------------------------------------------

int FieldDescr::insert_field(const std::string & field_name) throw()
{

  int id = field_count();

  // Check if field has already been inserted

  for (int i=0; i<id; i++) {
    if (field_name_[i] == field_name) {
      char buffer [ ERROR_LENGTH ];
      sprintf (buffer,
	       "Insert field called multiple times with same field %s",
	       field_name.c_str());
      WARNING("FieldDescr::insert_field", buffer);
      return i;
    }
  }
  // Insert field name and id

  field_name_.push_back(field_name);

  field_id_[field_name] = id;

  // Initialize attributes with default values

  precision_type precision = default_precision;

  bool * centered = new bool[3];
  centered[0] = true;
  centered[1] = true;
  centered[2] = true;

  int * ghosts = new int [3];
  ghosts[0] = 1;
  ghosts[1] = 1;
  ghosts[2] = 1;

  // int_set_type a;
  // field_in_group_.push_back(a);
  precision_. push_back(precision);
  centering_. push_back(centered);
  ghosts_.    push_back(ghosts);
  min_value_. push_back(0);
  max_value_. push_back(0);
  min_action_.push_back(field_action_none);
  max_action_.push_back(field_action_none);

  return id;
}

//----------------------------------------------------------------------

void FieldDescr::insert_group(const std::string & group_name) throw()
{
  int id = group_count();
  group_name_.push_back(group_name);
  group_id_[group_name] = id;
}

//----------------------------------------------------------------------

// void FieldDescr::set_field_in_group(int id_field, int id_group) 
//   throw(std::out_of_range)
// {
//   field_in_group_.at(id_field).insert(id_group);
// }

//----------------------------------------------------------------------

void FieldDescr::set_alignment(int alignment) throw()
{
  alignment_ = alignment;
}

//----------------------------------------------------------------------

void FieldDescr::set_padding(int padding) throw()
{
  padding_ = padding;
}

//----------------------------------------------------------------------

void FieldDescr::set_courant(double courant) throw()
{
  courant_ = courant;
}

//----------------------------------------------------------------------

void FieldDescr::set_precision(int id_field, precision_type precision) 
  throw(std::out_of_range)
{
  if ( ! cello::is_precision_supported (precision) ) {
    char buffer[80];
    sprintf (buffer,"precision \"%s\" is not supported",
	     cello::precision_name[precision]);
    WARNING("FieldDescr::set_precision",
	    buffer);
  }
  precision_.at(id_field) = 
    (precision == precision_default) ? default_precision : precision;
}

//----------------------------------------------------------------------

int FieldDescr::bytes_per_element(int id_field) const throw()
{
  return cello::sizeof_precision (precision(id_field));
}

//----------------------------------------------------------------------

void FieldDescr::set_centering(int id_field, bool cx, bool cy, bool cz) 
  throw(std::out_of_range)
{
  centering_.at(id_field)[0] = cx;
  centering_.at(id_field)[1] = cy;
  centering_.at(id_field)[2] = cz;
}

//----------------------------------------------------------------------

void FieldDescr::set_ghosts(int id_field, int gx, int gy, int gz) 
  throw(std::out_of_range)
{
  ghosts_.at(id_field)[0] = gx;
  ghosts_.at(id_field)[1] = gy;
  ghosts_.at(id_field)[2] = gz;
}

//----------------------------------------------------------------------

void FieldDescr::set_minimum 
(
 int               id_field, 
 double            min_value, 
 field_action_type min_action
 ) throw(std::out_of_range)
{
  min_value_.at(id_field)  = min_value;
  min_action_.at(id_field) = min_action;
}

//----------------------------------------------------------------------

void FieldDescr::set_maximum 
(
 int               id_field, 
 double            max_value, 
 field_action_type max_action
 ) throw(std::out_of_range)
{
  max_value_.at(id_field)  = max_value;
  max_action_.at(id_field) = max_action;
}

//======================================================================

void FieldDescr::copy_(const FieldDescr & field_descr) throw()
{
  alignment_      = field_descr.alignment_;
  padding_        = field_descr.padding_;
  courant_        = field_descr.courant_;
  field_name_     = field_descr.field_name_;
  field_id_       = field_descr.field_id_;
  group_name_     = field_descr.group_name_;
  group_id_       = field_descr.group_id_;
  // field_in_group_ = field_descr.field_in_group_;
  precision_      = field_descr.precision_;
  centering_.clear();
  for (size_t i=0; i<field_descr.centering_.size(); i++) {
    centering_.push_back(new bool[3]);
    centering_.at(i)[0] = field_descr.centering_.at(i)[0];
    centering_.at(i)[1] = field_descr.centering_.at(i)[1];
    centering_.at(i)[2] = field_descr.centering_.at(i)[2];
  }
  ghosts_.clear();
  for (size_t i=0; i<field_descr.ghosts_.size(); i++) {
    ghosts_.push_back(new int[3]);
    ghosts_.at(i)[0] = field_descr.ghosts_.at(i)[0];
    ghosts_.at(i)[1] = field_descr.ghosts_.at(i)[1];
    ghosts_.at(i)[2] = field_descr.ghosts_.at(i)[2];
  }
  min_value_      = field_descr.min_value_;
  max_value_      = field_descr.max_value_;
  min_action_     = field_descr.min_action_;
  max_action_     = field_descr.max_action_;
}

