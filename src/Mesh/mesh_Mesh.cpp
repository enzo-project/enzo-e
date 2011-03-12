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
 GroupProcess * group_process,
 int nx, int ny, int nz,
 int nbx, int nby, int nbz
 ) throw ()
  : group_process_(group_process),
    tree_(0),
    dimension_(0),
    refine_(0),
    max_level_(0),
    balanced_(0),
    backfill_(0),
    coalesce_(0)

{
  // Initialize extents
  for (int i=0; i<3; i++) {
    lower_[i] = 0.0;
    upper_[i] = 1.0;
  }

}

//----------------------------------------------------------------------

Mesh::~Mesh() throw()
{
  for (size_t i=0; i<patch_list_.size(); i++) {
    delete patch_list_[i];
    patch_list_[i] = 0;
  }
  delete tree_;
}

//----------------------------------------------------------------------

size_t Mesh::num_patches() const  throw()
{
  return patch_list_.size();
}

//----------------------------------------------------------------------

Patch * Mesh::patch(size_t i) const throw()
{
  return ( patch_list_.size()-1 >= i )? patch_list_[i] : 0;
}

//----------------------------------------------------------------------
// virtual
void Mesh::insert_patch(Patch * patch) throw()
{
  int size = patch_list_.size();
  patch_list_.resize(size + 1);
  patch_list_[size] = patch;
}
