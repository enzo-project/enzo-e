// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MappingArray.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-22
/// @brief    Mapping of Charm++ array Index to processors

#include "charm.hpp"

//======================================================================

MappingArray::MappingArray(int nx, int ny, int nz)
  :  CkArrayMap()
{
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
}

//----------------------------------------------------------------------

int MappingArray::procNum(int, const CkArrayIndex &idx) {

  int v3[3];

  v3[0] = idx.data()[0];
  v3[1] = idx.data()[1];
  v3[2] = idx.data()[2];

  Index in;
  in.set_values(v3);

  int ix,iy,iz;
  in.array    (&ix,&iy,&iz);

  int index = (ix + nx_*(iy + ny_*iz)) % CkNumPes();

  return index;
}

