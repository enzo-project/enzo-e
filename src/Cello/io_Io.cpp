// See LICENSE_CELLO file for license and copyright information

/// @file     io_Io.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

Io::Io() throw()
  : meta_name_()
{}

//----------------------------------------------------------------------

void Io::meta_value 
(int index, 
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{
  ASSERT1 ("Io::meta_value()",
   	   "index %d out of range", index,
   	   0 <= index && (size_t)index < meta_name_.size());

  if (name) (*name) = meta_name_[index];
  if (nxd) (*nxd) = 1;
  if (nyd) (*nyd) = 0;
  if (nzd) (*nzd) = 0;
}

//----------------------------------------------------------------------
void Io::field_array 
(void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
}

//----------------------------------------------------------------------

void Io::particle_array 
(int it, int ib, int ia,
 void ** buffer, std::string * name, int * type,
 int * n, int * k) throw()
{
}

//======================================================================

int Io::data_size () const
{
  int size = 0;

  int list_size = meta_name_.size();
  SIZE_SCALAR_TYPE(size,int,list_size);
  for (int i=0; i<list_size; i++) {
    SIZE_STRING_TYPE(size,meta_name_[i]);
  }

  return size;
}

//----------------------------------------------------------------------

char * Io::save_data (char * buffer) const
{
  char * pc = buffer;

  int list_size = meta_name_.size();
  SAVE_SCALAR_TYPE(pc,int,list_size);
  for (int i=0; i<list_size; i++) {
    SAVE_STRING_TYPE(pc,meta_name_[i]);
  }

  ASSERT2 ("Io::save_data()",
  	   "Expecting buffer size %d actual size %d",
  	   Io::data_size(),(pc-buffer),
  	   (Io::data_size() == (pc-buffer)));
  
  // return first byte after filled buffer
  return pc;
}

//----------------------------------------------------------------------

char * Io::load_data (char * buffer)
{
  char * pc = buffer;

  int list_size = 0;
  LOAD_SCALAR_TYPE(pc,int,list_size);
  meta_name_.clear();
  meta_name_.resize(list_size);
  for (int i=0; i<list_size; i++) {
    LOAD_STRING_TYPE(pc,meta_name_[i]);
  }

  return pc;
}

