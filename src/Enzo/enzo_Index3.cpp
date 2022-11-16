// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_Index3.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-10-05
///
/// @brief Implementation of Index3 class for indexing a sparse 3D array

#include "enzo.hpp"

// #define DEBUG_INDEX

//----------------------------------------------------------------------

Index3::Index3()
  : v_ {0,0,0}
{
}

//----------------------------------------------------------------------

Index3::Index3(int ix, int iy, int iz) 
{
  set_values(ix,iy,iz);
}

//----------------------------------------------------------------------

bool Index3::operator == (const Index3 & index) const
{
  return (v_[0] == index.v_[0] && 
	  v_[1] == index.v_[1] &&
	  v_[2] == index.v_[2]);
}

//----------------------------------------------------------------------

bool Index3::operator != (const Index3 & index) const
{
  return ! (*this == index);
}

// ----------------------------------------------------------------------

int Index3::data_size () const
{
  int size = 0;
  SIZE_ARRAY_TYPE(size,int,v_,3);
  return size;
}

// ----------------------------------------------------------------------

char * Index3::save_data (char * buffer) const
{
  char * pc = buffer;
  SAVE_ARRAY_TYPE(pc,int,v_,3);
  return pc;
}

// ----------------------------------------------------------------------

char * Index3::load_data (char * buffer)
{
  char * pc = buffer;
  LOAD_ARRAY_TYPE(pc,int,v_,3);
  return pc;
}

