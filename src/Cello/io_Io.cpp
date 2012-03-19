// See LICENSE_CELLO file for license and copyright information

/// @file     io_Io.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

Io::Io(size_t meta_count, size_t data_count) throw()
  : meta_count_(meta_count), 
    meta_name_(), 
    data_count_(data_count),
    data_name_()
{}

//----------------------------------------------------------------------

void Io::meta_value 
(int index, 
 void ** buffer, std::string * name, enum scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
  ASSERT1 ("Io::meta_value()",
   	   "index %d out of range", index,
   	   0 <= index && (size_t)index < meta_count_);

  if (name) (*name) = meta_name_[index];
  if (nxd) (*nxd) = 1;
  if (nyd) (*nyd) = 0;
  if (nzd) (*nzd) = 0;
}

//----------------------------------------------------------------------
void Io::data_value 
(int index, 
 void ** buffer, std::string * name, enum scalar_type * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
}

//======================================================================

