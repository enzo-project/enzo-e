// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoHierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoHierarchy class

#include "io.hpp"

//----------------------------------------------------------------------

IoHierarchy::IoHierarchy(const Hierarchy * hierarchy) throw ()
  : Io(),
    hierarchy_((Hierarchy *)hierarchy)
{
  meta_name_.push_back("lower");
  meta_name_.push_back("upper");
  meta_name_.push_back("max_level");
}


//----------------------------------------------------------------------

void IoHierarchy::meta_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
  Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);

  int count = 0;

  if (index == count++) {

    *buffer = (void *) hierarchy_->lower_;
    *type   = scalar_type_double;
    *nxd     = 3;
    
  } else if (index == count++) {

    *buffer = (void *) hierarchy_->upper_;
    *type   = scalar_type_double;
    *nxd     = 3;

  } else if (index == count++) {

    *buffer = (void *) & hierarchy_->max_level_;
    *type   = scalar_type_int;

  }
}

//----------------------------------------------------------------------

void IoHierarchy::data_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
}

//----------------------------------------------------------------------
