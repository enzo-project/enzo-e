// See LICENSE_CELLO file for license and copyright information

/// @file     Io_IoHierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoHierarchy class

#include "io.hpp"

//----------------------------------------------------------------------

IoHierarchy::IoHierarchy(const Hierarchy * hierarchy) throw ()
  : Io(3,0),
    hierarchy_(hierarchy)
{
  meta_name_.push_back("lower");
  meta_name_.push_back("upper");
  meta_name_.push_back("patch_count");
}


//----------------------------------------------------------------------

void IoHierarchy::meta_value
(int index,
 void ** buffer, const char ** name, enum scalar_type * type,
 int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  Io::meta_value(index,buffer,name,type,n0,n1,n2,n3,n4);

  if (index == 0) {

    *buffer = (void *) hierarchy_->lower_;
    *type   = scalar_type_double;
    *n0     = 3;
    
  } else if (index == 1) {

    *buffer = (void *) hierarchy_->upper_;
    *type   = scalar_type_double;
    *n0     = 3;

  } else if (index == 2) {

    *buffer = (void *) & hierarchy_->patch_count_;
    *type   = scalar_type_int;
  }
  
}

//----------------------------------------------------------------------

void IoHierarchy::data_value
(int index,
 void ** buffer, const char ** name, enum scalar_type * type,
 int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
}

//----------------------------------------------------------------------
