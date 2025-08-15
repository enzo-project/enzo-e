// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-14
/// @brief    Implementation of the FieldDescr class

#include "cello.hpp"

#include "data.hpp"

//----------------------------------------------------------------------

FieldDescr::FieldDescr () throw ()
  : name_(),
    num_permanent_(0),
    num_temporary_(0),
    id_(),
    groups_(),
    alignment_(1),
    padding_(0),
    precision_(),
    centering_(),
    ghost_depth_(),
    conserved_(),
    history_(0),
    history_id_()
{
  for (int i=0; i<3; i++) {
    ghost_depth_default_[i] = 0;
  }
}

//----------------------------------------------------------------------

FieldDescr::~FieldDescr() throw()
{
}

//----------------------------------------------------------------------

FieldDescr::FieldDescr(const FieldDescr & field_descr) throw() 
{
  copy_(field_descr);
}

//----------------------------------------------------------------------

FieldDescr & FieldDescr::operator= (const FieldDescr & field_descr) throw()
{
  copy_(field_descr);

  return *this;
}

//----------------------------------------------------------------------

int FieldDescr::field_count() const throw()
{
  return name_.size();
}

//----------------------------------------------------------------------

std::string FieldDescr::field_name(int id_field) const throw()
{
  ASSERT1("FieldDescr::field_name",
          "index %d is out of range for vector of field names.", id_field,
          name_.size() > (std::vector<std::string>::size_type)id_field);
  return name_[id_field];
}

//----------------------------------------------------------------------

bool FieldDescr::is_field(const std::string & name) const throw()
{ 
  return (id_.find(name) != id_.end());
}

//----------------------------------------------------------------------

int FieldDescr::field_id(const std::string & name) const throw()
{
  auto it = id_.find(name);
  if (it != id_.end()) {
    return it->second;
  } else {
    //    WARNING1("FieldDescr::field_id()",
    //	   "Trying to access unknown Field \"%s\"",
    //	   name.c_str());
    return -1;
  }
}

//----------------------------------------------------------------------

void FieldDescr::centering
(
 int id_field,
 int * cx, 
 int * cy, 
 int * cz
 ) const throw()
{
  if (id_field>=0) {
    if (cx) (*cx) = std::get<0>(centering_[id_field]);
    if (cy) (*cy) = std::get<1>(centering_[id_field]);
    if (cz) (*cz) = std::get<2>(centering_[id_field]);
  }
}

//----------------------------------------------------------------------

void FieldDescr::ghost_depth
(
 int id_field,
 int * gx, 
 int * gy, 
 int * gz
 ) const throw()
{
  if (id_field>=0) {
    int g3[3] = {
      std::get<0>(ghost_depth_[id_field]),
      std::get<1>(ghost_depth_[id_field]),
      std::get<2>(ghost_depth_[id_field])
    };
    int gd[3] = {
      ghost_depth_default_[0],
      ghost_depth_default_[1],
      ghost_depth_default_[2]
    };
    if (gx) (*gx) = g3[0] < 0 ? gd[0] : g3[0];
    if (gy) (*gy) = g3[1] < 0 ? gd[1] : g3[1];
    if (gz) (*gz) = g3[2] < 0 ? gd[2] : g3[2];
  }
}

//----------------------------------------------------------------------

int FieldDescr::insert_permanent(const std::string & field_name) throw()
{

  bool permanent;

  int id = insert_(field_name, permanent = true);

  ++ num_permanent_;

  return id;
}

//----------------------------------------------------------------------

int FieldDescr::insert_temporary(const std::string & field_name) throw()
{
  bool permanent;
  
  int id = insert_(field_name, permanent = false);

  ++ num_temporary_;

  return id;
}

//----------------------------------------------------------------------

int FieldDescr::insert_(const std::string & field_name,
			bool is_permanent) throw()
{
  // Assumes all permanent added before any temporary
  
  int id = num_permanent_ + num_temporary_;

  // Check if field has already been inserted

  if (is_permanent) {
    for (int i=0; i<id; i++) {
      if (name_[i] == field_name) {
	char buffer [ ERROR_LENGTH ];
	snprintf (buffer, ERROR_LENGTH,
		 "Insert field called multiple times with same field %s",
		 field_name.c_str());
	WARNING("FieldDescr::insert_permanent", buffer);
	return i;
      }
    }
  }

  // Save field name (unless anonymous temporary)
  
  if (field_name != "") {
    name_.push_back(field_name);
    id_[field_name] = id;
  }
  
  // Initialize attributes with default values

  int precision = default_precision;

  precision_.  push_back(precision);
  centering_.  push_back(std::tuple<int,int,int>(0,0,0));
  ghost_depth_.push_back(std::tuple<int,int,int>(-1,-1,-1));

  return id;
}

//----------------------------------------------------------------------

void FieldDescr::set_precision(int id_field, int precision) throw()
{
  if ( ! cello::is_precision_supported (precision) ) {
    WARNING1("FieldDescr::set_precision","precision \"%s\" is not supported",
	         cello::precision_name[precision]);
  }
  if (id_field >= 0) {
    precision_.at(id_field) = 
      (precision == precision_default) ? default_precision : precision;
  }
}

//----------------------------------------------------------------------

int FieldDescr::bytes_per_element(int id_field) const throw()
{
  return cello::sizeof_precision (precision(id_field));
}

//----------------------------------------------------------------------

void FieldDescr::set_centering(int id_field, int cx, int cy, int cz) throw()
{
  if (id_field >= 0) {
    std::get<0>(centering_[id_field]) = cx;
    std::get<1>(centering_[id_field]) = cy;
    std::get<2>(centering_[id_field]) = cz;
  }
}

//----------------------------------------------------------------------

void FieldDescr::set_ghost_depth(int id_field, int gx, int gy, int gz) throw()
{
  if (id_field >= 0) {
    std::get<0>(ghost_depth_[id_field]) = gx;
    std::get<1>(ghost_depth_[id_field]) = gy;
    std::get<2>(ghost_depth_[id_field]) = gz;
  }
}

//----------------------------------------------------------------------

void FieldDescr::set_default_ghost_depth(int gx, int gy, int gz) throw()
{
  ghost_depth_default_[0] = gx;
  ghost_depth_default_[1] = gy;
  ghost_depth_default_[2] = gz;
}

//======================================================================

void FieldDescr::copy_(const FieldDescr & field_descr) throw()
{
  name_      = field_descr.name_;
  num_permanent_ = field_descr.num_permanent_;
  num_temporary_ = field_descr.num_temporary_;
  id_        = field_descr.id_;
  groups_    = field_descr.groups_;
  alignment_ = field_descr.alignment_;
  padding_   = field_descr.padding_;
  precision_ = field_descr.precision_;
  centering_ = field_descr.centering_;
  ghost_depth_ = field_descr.ghost_depth_;
  for (int i=0; i<3; i++) {
    ghost_depth_default_[i] = field_descr.ghost_depth_default_[i];
  }
  conserved_.resize(field_descr.conserved_.size());
  for (size_t i=0; i<field_descr.conserved_.size(); i++) {
    conserved_[i] = field_descr.conserved_[i];
  }
  history_ = field_descr.history_;
  history_id_.resize(field_descr.history_id_.size());
  for (size_t i=0; i<field_descr.history_id_.size(); i++) {
    history_id_[i] = field_descr.history_id_[i];
  }
}
