// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-14
/// @brief    Implementation of the FieldDescr class

#include "cello.hpp"

#include "data.hpp"

//----------------------------------------------------------------------

FieldDescr::FieldDescr () 
  throw ()
  : name_(),
    num_permanent_(0),
    id_(),
    groups_(),
    alignment_(1),
    padding_(0),
    courant_(1),
    precision_(),
    centering_(),
    ghost_depth_()
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
  for (size_t i=0; i<ghost_depth_.size();    i++) {
    delete [] ghost_depth_[i];
    ghost_depth_[i] = 0;
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
  return name_.size();
}

//----------------------------------------------------------------------

std::string FieldDescr::field_name(size_t id_field) const 
  throw(std::out_of_range)
{ 
  return name_.at(id_field);
}

//----------------------------------------------------------------------

bool FieldDescr::is_field(const std::string & name) const 
  throw()
{ 
  return (id_.find(name) != id_.end());
}

//----------------------------------------------------------------------

int FieldDescr::field_id(const std::string & name) const
  throw()
{
  //  return id_[name]; // ERROR IN PGI ON GORDON 11.9-0 64-bit
  std::map<const std::string,int>::const_iterator it;
  it=id_.find(name);
  if (it != id_.end()) {
    return it->second;
  } else {
    WARNING1("FieldDescr::field_id()",
	   "Trying to access unknown Field \"%s\"",
	   name.c_str());
    return -1;
  }
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

double FieldDescr::courant() const throw()
{
  return courant_;
}

//----------------------------------------------------------------------

void FieldDescr::centering
(
 int id_field,
 int * cx, 
 int * cy, 
 int * cz
 ) const throw(std::out_of_range)
{
  if (cx) (*cx) = centering_.at(id_field)[0];
  if (cy) (*cy) = centering_.at(id_field)[1];
  if (cz) (*cz) = centering_.at(id_field)[2];
}

//----------------------------------------------------------------------

void FieldDescr::ghost_depth
(
 int id_field,
 int * gx, 
 int * gy, 
 int * gz
 ) const throw(std::out_of_range)
{
  if (gx) (*gx) = ghost_depth_.at(id_field)[0];
  if (gy) (*gy) = ghost_depth_.at(id_field)[1];
  if (gz) (*gz) = ghost_depth_.at(id_field)[2];
}

//----------------------------------------------------------------------

precision_type FieldDescr::precision(int id_field) const 
  throw(std::out_of_range)
{
  return precision_.at(id_field);
}

//----------------------------------------------------------------------

int FieldDescr::insert_permanent(const std::string & field_name) throw()
{

  int id = insert_(field_name);

  num_permanent_ = id+1;

  return id;
}

//----------------------------------------------------------------------

int FieldDescr::insert_temporary(const std::string & field_name) throw()
{
  return insert_(field_name);
}

//----------------------------------------------------------------------

int FieldDescr::insert_(const std::string & field_name) throw()
{
  int id = field_count();

  // Check if field has already been inserted

  for (int i=0; i<id; i++) {
    if (name_[i] == field_name) {
      char buffer [ ERROR_LENGTH ];
      sprintf (buffer,
	       "Insert field called multiple times with same field %s",
	       field_name.c_str());
      WARNING("FieldDescr::insert_permanent", buffer);
      return i;
    }
  }
  // Insert field name and id

  name_.push_back(field_name);

  id_[field_name] = id;

  // Initialize attributes with default values

  precision_type precision = default_precision;

  int * centered = new int[3];
  centered[0] = 0;
  centered[1] = 0;
  centered[2] = 0;

  int * ghost_depth = new int [3];
  ghost_depth[0] = 1;
  ghost_depth[1] = 1;
  ghost_depth[2] = 1;

  // int_set_type a;
  // field_in_group_.push_back(a);
  precision_. push_back(precision);
  centering_. push_back(centered);
  ghost_depth_.    push_back(ghost_depth);

  return id;
}

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

void FieldDescr::set_centering(int id_field, int cx, int cy, int cz) 
  throw(std::out_of_range)
{
  centering_.at(id_field)[0] = cx;
  centering_.at(id_field)[1] = cy;
  centering_.at(id_field)[2] = cz;
}

//----------------------------------------------------------------------

void FieldDescr::set_ghost_depth(int id_field, int gx, int gy, int gz) 
  throw(std::out_of_range)
{
  ghost_depth_.at(id_field)[0] = gx;
  ghost_depth_.at(id_field)[1] = gy;
  ghost_depth_.at(id_field)[2] = gz;
}

//======================================================================

void FieldDescr::copy_(const FieldDescr & field_descr) throw()
{
  alignment_ = field_descr.alignment_;
  padding_   = field_descr.padding_;
  courant_   = field_descr.courant_;
  name_      = field_descr.name_;
  num_permanent_ = field_descr.num_permanent_;
  id_        = field_descr.id_;
  groups_    = field_descr.groups_;
  precision_ = field_descr.precision_;
  centering_.clear();
  for (size_t i=0; i<field_descr.centering_.size(); i++) {
    centering_.push_back(new int[3]);
    centering_.at(i)[0] = field_descr.centering_.at(i)[0];
    centering_.at(i)[1] = field_descr.centering_.at(i)[1];
    centering_.at(i)[2] = field_descr.centering_.at(i)[2];
  }
  ghost_depth_.clear();
  for (size_t i=0; i<field_descr.ghost_depth_.size(); i++) {
    ghost_depth_.push_back(new int[3]);
    ghost_depth_.at(i)[0] = field_descr.ghost_depth_.at(i)[0];
    ghost_depth_.at(i)[1] = field_descr.ghost_depth_.at(i)[1];
    ghost_depth_.at(i)[2] = field_descr.ghost_depth_.at(i)[2];
  }
}


