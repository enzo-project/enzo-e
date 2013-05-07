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

  int vx = idx.data()[0];
  int vy = idx.data()[1];
  int vz = idx.data()[2];

  Index in;
  int ix,iy,iz;
  in.set_values(vx, vy, vz);
  in.array    (&ix,&iy,&iz);

  TRACE3("ArrayMap  newindex = %d %d %d",ix,iy,iz);

  
  TRACE3("ArrayMap  size  = %d %d %d",nx_,ny_,nz_);

  int index = (ix + nx_*(iy + ny_*iz)) % CkNumPes();
  TRACE1("ArrayMap  proc = %d",index);

  return index;
}

