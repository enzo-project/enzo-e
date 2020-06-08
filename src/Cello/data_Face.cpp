// See LICENSE_CELLO file for license and copyright information

/// @file     data_Face.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-22
/// @brief    Implementation of the Face class

#include "data.hpp"

// #define DEBUG_FACE

// #define DEBUG_REFRESH

//======================================================================

int Face::data_size () const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH Face::data_size()\n",CkMyPe());
#endif  
  int size = 0;

  size += 3*sizeof(int);  // int ix_,iy_,iz_;
  size += 3*sizeof(int);  // int rx_, ry_, rz_;
}

//----------------------------------------------------------------------

char * Face::save_data (char * buffer) const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH Face::save_data()\n",CkMyPe());
#endif  
  union {
    char  * pc;
    int   * pi;
    float * pf;
  };

  pc = (char *) buffer;

  (*pi++) = ix_;
  (*pi++) = iy_;
  (*pi++) = iz_;

  (*pi++) = rx_;
  (*pi++) = ry_;
  (*pi++) = rz_;
  
  ASSERT2("Face::save_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------

char * Face::load_data (char * buffer)
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH Face::load_data()\n",CkMyPe());
#endif  
  union {
    int   * pi;
    char  * pc;
    float * pf;
  };

  pc = (char *) buffer;

  ix_ = (*pi++);
  iy_ = (*pi++);
  iz_ = (*pi++);

  rx_ = (*pi++);
  ry_ = (*pi++);
  rz_ = (*pi++);

  ASSERT2("Face::load_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

