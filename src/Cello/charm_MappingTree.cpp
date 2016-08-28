// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MappingTree.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-22
/// @brief    Mapping of Charm++ array Index to processors

#include "charm.hpp"

//======================================================================

MappingTree::MappingTree(int nx, int ny, int nz)
  :  CkArrayMap()
{
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
}

//----------------------------------------------------------------------

int MappingTree::procNum(int, const CkArrayIndex &idx) {

  int v3[3];

  v3[0] = idx.data()[0];
  v3[1] = idx.data()[1];
  v3[2] = idx.data()[2];

  Index in;
  in.set_values(v3);

  const int level = in.level();

  int iax,iay,iaz;
  int itx,ity,itz;
  in.array (&iax,&iay,&iaz);
  in.tree  (&itx,&ity,&itz);
  
  int ix=iax;
  int iy=iay;
  int iz=iaz;
  for (int l=0; l<level; l++) {
    int icx=0,icy=0,icz=0;
    in.child(l+1,&icx,&icy,&icz);
    ix = (ix<<1) | icx;
    iy = (iy<<1) | icy;
    iz = (iz<<1) | icz;
  }

  int index = (ix + nx_*(iy + ny_*iz)) % CkNumPes();

  return index;
}

