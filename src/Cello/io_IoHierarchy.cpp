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
  meta_name_.push_back("lower");
  meta_name_.push_back("upper");
  meta_name_.push_back("max_level");
  int i;
  for (i=0; i<3; i++) {
    lower_[i] = hierarchy->lower_[i];
    upper_[i] = hierarchy->upper_[i];
  }
  max_level_ = hierarchy->max_level_;
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

    *buffer = (void *) lower_;
    *type   = type_double;
    *nxd     = 3;
    
  } else if (index == count++) {

    *buffer = (void *) upper_;
    *type   = type_double;
    *nxd     = 3;

  } else if (index == count++) {

    *buffer = (void *) & max_level_;
    *type   = type_int;

  }
}

//----------------------------------------------------------------------

void IoHierarchy::save_to (void * v)
{
  Hierarchy * h = static_cast<Hierarchy*> (v);

  for (int i=0; i<3; i++) {
    h->lower_[i] = lower_[i];
    h->upper_[i] = upper_[i];
  }
  h->max_level_ = max_level_;
}
