// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoHierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoHierarchy class

#include "io.hpp"

//----------------------------------------------------------------------

IoHierarchy::IoHierarchy(const Hierarchy * hierarchy) throw ()
  : Io()
{
  meta_name_.push_back("max_level");
  meta_name_.push_back("lower");
  meta_name_.push_back("upper");
  meta_name_.push_back("root_size");
  meta_name_.push_back("blocking");

  max_level_ = hierarchy->max_level_;
  int i;
  for (i=0; i<3; i++) {
    lower_[i] = hierarchy->lower_[i];
    upper_[i] = hierarchy->upper_[i];
    root_size_[i] = hierarchy->root_size_[i];
    blocking_[i] = hierarchy->blocking_[i];
  }
}


//----------------------------------------------------------------------

void IoHierarchy::meta_value
(int index,
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{
  Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);

  int count = 0;

  if (index == count++) {

    *buffer = (void *) & max_level_;
    *type   = type_int;

  } else if (index == count++) {

    *buffer = (void *) lower_;
    *type   = type_double;
    *nxd     = 3;

  } else if (index == count++) {

    *buffer = (void *) upper_;
    *type   = type_double;
    *nxd     = 3;

  } else if (index == count++) {

    *buffer = (void *) root_size_;
    *type   = type_int;
    *nxd    = 3;

  } else if (index == count++) {

    *buffer = (void *) blocking_;
    *type   = type_int;
    *nxd    = 3;

  }
}

//----------------------------------------------------------------------

void IoHierarchy::save_to (void * v)
{
  Hierarchy * h = static_cast<Hierarchy*> (v);

  h->max_level_ = max_level_;
  for (int i=0; i<3; i++) {
    h->lower_[i] = lower_[i];
    h->upper_[i] = upper_[i];
    h->root_size_[i] = root_size_[i];
    h->blocking_[i] = blocking_[i];
  }
}
