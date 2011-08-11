// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Brief description of file mesh_Mesh.cpp

#include "cello.hpp"

#include "mesh.hpp"

#include "mesh_charm.hpp"

//----------------------------------------------------------------------

Mesh::Mesh
(
 const Factory * factory,
 int nx, int ny, int nz,
 int nbx, int nby, int nbz
 ) throw ()
  : factory_(factory)
    // tree_(0),
    // dimension_(0),
    // refine_(0),
    // max_level_(0),
    // balanced_(0),
    // backfill_(0),
    // coalesce_(0)

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

void Mesh::set_lower(double x, double y, double z) throw ()
{
  lower_[0] = x;
  lower_[1] = y;
  lower_[2] = z;
}

//----------------------------------------------------------------------

void Mesh::set_upper(double x, double y, double z) throw ()
{
  upper_[0] = x;
  upper_[1] = y;
  upper_[2] = z;
}

//----------------------------------------------------------------------

// int Mesh::dimension() const throw ()
// { 
//   return dimension_; 
// }

// //----------------------------------------------------------------------

// void Mesh::set_dimension(int dimension) throw ()
// {
//   dimension_ = dimension; 
// }

//----------------------------------------------------------------------

int Mesh::dimension() const throw ()
{
  if (patch_list_.size() == 0) {
    return 0;
  } else {
    int nx,ny,nz;
    patch_list_[0]->size(&nx,&ny,&nz);
    if (nz != 1) return 3;
    if (ny != 1) return 2;
    if (nx != 1) return 1;
    return 0;
  }
}

//----------------------------------------------------------------------

void Mesh::lower(double * x, double * y, double * z) const throw ()
{
  if (x) *x = lower_[0];
  if (y) *y = lower_[1];
  if (z) *z = lower_[2];
}

//----------------------------------------------------------------------

void Mesh::upper(double * x, double * y, double * z) const throw ()
{
  if (x) *x = upper_[0];
  if (y) *y = upper_[1];
  if (z) *z = upper_[2];
}
// //----------------------------------------------------------------------

// int Mesh::max_level() const throw ()
// { 
//   return max_level_; 
// }

// //----------------------------------------------------------------------

// void Mesh::set_max_level(int max_level) throw ()
// {
//   max_level_ = max_level; 
// }

// //----------------------------------------------------------------------

// int Mesh::refine_factor() const throw ()
// {
//   return refine_; 
// }

// //----------------------------------------------------------------------

// void Mesh::set_refine_factor(int refine) throw ()
// {
//   refine_ = refine; 
// }

//----------------------------------------------------------------------

// Patch * Mesh::root_patch() throw ()
// { 
//   return (patch_list_.size() > 0) ? patch_list_[0] : NULL; 
// }

//----------------------------------------------------------------------

size_t Mesh::num_patches() const throw()
{
  return patch_list_.size();
}

//----------------------------------------------------------------------

Patch * Mesh::patch(size_t i) throw()
{
  return ( patch_list_.size()-1 >= i )? patch_list_[i] : 0;
}

//----------------------------------------------------------------------

Patch * Mesh::patch(size_t i) const throw()
{
  return ( patch_list_.size()-1 >= i )? patch_list_[i] : 0;
}

//----------------------------------------------------------------------

void Mesh::insert_patch(Patch * patch) throw()
{
  int size = patch_list_.size();
  patch_list_.resize(size + 1);
  patch_list_[size] = patch;
}

// //----------------------------------------------------------------------

// bool Mesh::balanced() const throw ()
// {
//   return balanced_; 
// }

// //----------------------------------------------------------------------

// void Mesh::set_balanced(bool balanced) throw ()
// { 
//   balanced_ = balanced; 
// }

// //----------------------------------------------------------------------

// bool Mesh::backfill() const throw ()
// {
//   return backfill_; 
// }

// //----------------------------------------------------------------------

// void Mesh::set_backfill(bool backfill) throw ()
// { 
//   backfill_ = backfill; 
// }

// //----------------------------------------------------------------------

// bool Mesh::coalesce() const throw ()
// {
//   return coalesce_; 
// }

// //----------------------------------------------------------------------

// void Mesh::set_coalesce(bool coalesce) throw ()
// { 
//   coalesce_ = coalesce; 
// }

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
#  include "mesh.def.h"
#endif
