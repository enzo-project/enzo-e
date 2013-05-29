// See LICENSE_CELLO file for license and copyright information

/// @file     charm_ArrayMap.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-22
/// @brief    Mapping of Charm++ array Index to processors

#include "charm.hpp"

//======================================================================

ArrayMap::ArrayMap(int nx, int ny, int nz) 
  :  CBase_ArrayMap()
{
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
}

//----------------------------------------------------------------------

ArrayMap::~ArrayMap() {
}

int ArrayMap::procNum(int, const CkArrayIndex &idx) {

  int v3[3];

  v3[0] = idx.data()[0];
  v3[1] = idx.data()[1];
  v3[2] = idx.data()[2];

  Index in;
  int ix,iy,iz;
  in.set_values(v3);
  in.array    (&ix,&iy,&iz);

  TRACE3("ArrayMap  newindex = %d %d %d",ix,iy,iz);

  
  TRACE3("ArrayMap  size  = %d %d %d",nx_,ny_,nz_);

  int index = (ix + nx_*(iy + ny_*iz)) % CkNumPes();
  TRACE1("ArrayMap  proc = %d",index);

  return index;
}

