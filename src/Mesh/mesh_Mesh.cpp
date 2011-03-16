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
 Factory * factory,
 GroupProcess * group_process,
 int nx, int ny, int nz,
 int nbx, int nby, int nbz
 ) throw ()
  : factory_(factory),
    group_process_(group_process),
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

int Mesh::dimension() const throw ()
{ 
  return dimension_; 
}

//----------------------------------------------------------------------

void Mesh::set_dimension(int dimension) throw ()
{
  dimension_ = dimension; 
}

//----------------------------------------------------------------------

void Mesh::lower(double * nx, double * ny, double * nz) const throw ()
{
  *nx = lower_[0];
  *ny = lower_[1];
  *nz = lower_[2];
}

//----------------------------------------------------------------------

void Mesh::set_lower(double nx, double ny, double nz) throw ()
{
  lower_[0] = nx;
  lower_[1] = ny;
  lower_[2] = nz;
}

//----------------------------------------------------------------------

void Mesh::upper(double * nx, double * ny, double * nz) const throw ()
{
  *nx = upper_[0];
  *ny = upper_[1];
  *nz = upper_[2];
}

//----------------------------------------------------------------------

void Mesh::set_upper(double nx, double ny, double nz) throw ()
{
  upper_[0] = nx;
  upper_[1] = ny;
  upper_[2] = nz;
}

//----------------------------------------------------------------------

int Mesh::max_level() const throw ()
{ 
  return max_level_; 
}

//----------------------------------------------------------------------

void Mesh::set_max_level(int max_level) throw ()
{
  max_level_ = max_level; 
}

//----------------------------------------------------------------------

int Mesh::refine_factor() const throw ()
{
  return refine_; 
}

//----------------------------------------------------------------------

void Mesh::set_refine_factor(int refine) throw ()
{
  refine_ = refine; 
}

//----------------------------------------------------------------------

void Mesh::root_size(int * nx, int * ny, int * nz) const throw ()
{
  *nx = root_size_[0];
  *ny = root_size_[1];
  *nz = root_size_[2];
}

//----------------------------------------------------------------------

Patch * Mesh::root_patch() throw ()
{ 
  return (patch_list_.size() > 0) ? patch_list_[0] : 0; 
}

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

//----------------------------------------------------------------------

bool Mesh::balanced() const throw ()
{
  return balanced_; 
}

//----------------------------------------------------------------------

void Mesh::set_balanced(bool balanced) throw ()
{ 
  balanced_ = balanced; 
}

//----------------------------------------------------------------------

bool Mesh::backfill() const throw ()
{
  return backfill_; 
}

//----------------------------------------------------------------------

void Mesh::set_backfill(bool backfill) throw ()
{ 
  backfill_ = backfill; 
}

//----------------------------------------------------------------------

bool Mesh::coalesce() const throw ()
{
  return coalesce_; 
}

//----------------------------------------------------------------------

void Mesh::set_coalesce(bool coalesce) throw ()
{ 
  coalesce_ = coalesce; 
}

//----------------------------------------------------------------------

GroupProcess * Mesh::group()  const throw()
{
  return group_process_;
}


