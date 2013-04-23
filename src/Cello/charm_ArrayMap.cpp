// See LICENSE_CELLO file for license and copyright information

/// @file     charm_ArrayMap.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-22
/// @brief    Mapping of Charm++ array Index to processors

#include "charm.hpp"

//======================================================================

ArrayMap::ArrayMap(int nx, int ny, int nz) {
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
}

//----------------------------------------------------------------------

ArrayMap::~ArrayMap() {
}

int ArrayMap::procNum(int, const CkArrayIndex &idx) {

  const int ix = idx.data()[0];
  const int iy = idx.data()[1];
  const int iz = idx.data()[2];

  TRACE3("ArrayMap  index = %d %d %d",ix,iy,iz);
  TRACE3("ArrayMap  array = %d %d %d",nx_,ny_,nz_);
  int index = (ix + nx_*(iy + ny_*iz)) % CkNumPes();
  return index;
}

