// See LICENSE_CELLO file for license and copyright information

/// @file     io_Io.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

Io::Io(int meta_count, int data_count) throw()
    : meta_count_(meta_count), 
      meta_name_(), 
      data_count_(data_count),
      data_name_()
  {}

//----------------------------------------------------------------------

void Io::meta_value 
  (int index, 
   void ** buffer, const char ** name, enum scalar_type * type,
   int * n0, int * n1, int * n2, int * n3, int * n4) throw()
  {
    ASSERT1 ("Io::meta_value()",
	     "index %d out of range", index,
	     0 <= index && index < meta_count_);

    if (name) (*name) = meta_name_[index].c_str();
    if (n0) (*n0) = 1;
    if (n1) (*n1) = 0;
    if (n2) (*n2) = 0;
    if (n3) (*n3) = 0;
    if (n4) (*n4) = 0;
  }

//======================================================================

