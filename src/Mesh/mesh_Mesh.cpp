// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Brief description of file mesh_Mesh.cpp

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Mesh::Mesh
(
 DataDescr * data_descr,
 int nx, int ny, int nz,
 int nbx, int nby, int nbz
 ) throw ()
  : root_patch_(new Patch(nx,ny,nz,nbx,nby,nbz)),
    tree_(0),
    dimension_(0),
    refine_(0),
    max_level_(0),
    balanced_(0),
    backfill_(0),
    coalesce_(0)

{
  for (int i=0; i<3; i++) {
    lower_[i] = 0.0;
    upper_[i] = 1.0;
  }
  root_size_[0] = nx;
  root_size_[1] = ny;
  root_size_[2] = nz;
}

//----------------------------------------------------------------------

Mesh::~Mesh() throw()
{
  delete root_patch_;
  delete tree_;
}

