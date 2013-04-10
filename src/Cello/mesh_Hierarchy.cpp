// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Hierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Implementation of the Hierarchy class

#include "cello.hpp"

#include "mesh.hpp"

#include "charm_mesh.hpp"

//----------------------------------------------------------------------

Hierarchy::Hierarchy ( const Factory * factory,
		       int dimension, int refinement) throw ()
  : factory_((Factory * )factory),
    dimension_(dimension),
    refinement_(refinement),
    patch_count_(0),
    patch_tree_(new Tree (dimension,refinement))
{
  // Initialize extents
  for (int i=0; i<3; i++) {
    lower_[i] = 0.0;
    upper_[i] = 1.0;
    root_size_[i] = 1;
  }
}

//----------------------------------------------------------------------

Hierarchy::~Hierarchy() throw()
{
  ItNode it_node (patch_tree_);
  while (Node * node = it_node.next_leaf()) {
    Patch * patch = (Patch *) node->data();
    delete patch;
    patch = 0;
  }
  delete patch_tree_;
  patch_count_ = 0;
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
void Hierarchy::pup (PUP::er &p)
{
    
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  bool up = p.isUnpacking();

  p | factory_; // PUP::able
  p | dimension_;
  p | refinement_;
  p | patch_count_;
  if (up) patch_tree_ = new Tree (dimension_, refinement_);
  p | *patch_tree_;
  PUParray(p,root_size_,3);
  PUParray(p,lower_,3);
  PUParray(p,upper_,3);

}
#endif

//----------------------------------------------------------------------

void Hierarchy::set_lower(double x, double y, double z) throw ()
{
  lower_[0] = x;
  lower_[1] = y;
  lower_[2] = z;
}

//----------------------------------------------------------------------

void Hierarchy::set_upper(double x, double y, double z) throw ()
{
  upper_[0] = x;
  upper_[1] = y;
  upper_[2] = z;
}

//----------------------------------------------------------------------

void Hierarchy::set_root_size(int nx, int ny, int nz) throw ()
{
  root_size_[0] = nx;
  root_size_[1] = ny;
  root_size_[2] = nz;
}

//----------------------------------------------------------------------

// int Hierarchy::dimension() const throw ()
// { 
//   return dimension_; 
// }

// //----------------------------------------------------------------------

// void Hierarchy::set_dimension(int dimension) throw ()
// {
//   dimension_ = dimension; 
// }

//----------------------------------------------------------------------

int Hierarchy::dimension() const throw ()
{
  Patch * root = (Patch * )patch_tree_->root_node();
  if (root == 0) {
    return 0;
  } else {
    int nx,ny,nz;
    root->size(&nx,&ny,&nz);
    if (nz != 1) return 3;
    if (ny != 1) return 2;
    if (nx != 1) return 1;
    return 0;
  }
}

//----------------------------------------------------------------------

void Hierarchy::root_size(int * nx, int * ny, int * nz) const throw ()
{
  if (nx) *nx = root_size_[0];
  if (ny) *ny = root_size_[1];
  if (nz) *nz = root_size_[2];
}

//----------------------------------------------------------------------

void Hierarchy::lower(double * x, double * y, double * z) const throw ()
{
  if (x) *x = lower_[0];
  if (y) *y = lower_[1];
  if (z) *z = lower_[2];
}
//----------------------------------------------------------------------

void Hierarchy::upper(double * x, double * y, double * z) const throw ()
{
  if (x) *x = upper_[0];
  if (y) *y = upper_[1];
  if (z) *z = upper_[2];
}
// //----------------------------------------------------------------------

// int Hierarchy::max_level() const throw ()
// { 
//   return max_level_; 
// }

// //----------------------------------------------------------------------

// void Hierarchy::set_max_level(int max_level) throw ()
// {
//   max_level_ = max_level; 
// }

// //----------------------------------------------------------------------

// int Hierarchy::refine_factor() const throw ()
// {
//   return refine_; 
// }

// //----------------------------------------------------------------------

// void Hierarchy::set_refine_factor(int refine) throw ()
// {
//   refine_ = refine; 
// }

//----------------------------------------------------------------------

// Patch * Hierarchy::root_patch() throw ()
// { 
//   return (patch_list_.size() > 0) ? patch_list_[0] : NULL; 
// }

//----------------------------------------------------------------------

size_t Hierarchy::num_patches() const throw()
{
  return patch_tree_->num_nodes();
}

//----------------------------------------------------------------------

Patch * Hierarchy::patch(size_t i) throw()
{
  return (Patch * ) patch_tree_->root_node()->data();
}

//----------------------------------------------------------------------

Patch * Hierarchy::patch(size_t i) const throw()
{
  return (Patch * ) patch_tree_->root_node()->data();
}

//----------------------------------------------------------------------

void Hierarchy::create_root_patch 
(
 FieldDescr   * field_descr,
 int nx, int ny, int nz,
 int nbx, int nby, int nbz,
 bool allocate_blocks,
 int process_first, int process_last_plus) throw()
{

  // Create new empty patch

#ifdef CONFIG_USE_CHARM
  CProxy_Patch * root_patch;
#else
  Patch * root_patch;
#endif

  root_patch = factory()->create_patch
    (
     field_descr,
     nx,ny,nz,    // size
     0,0,0,       // offset
     nbx,nby,nbz, // blocking
     lower_[0], lower_[1], lower_[2],
     upper_[0], upper_[1], upper_[2],
     allocate_blocks,
     process_first, process_last_plus);

  patch_tree_->root_node()->set_data(root_patch);
  TRACE1 ("root_patch = %p",root_patch);
  {  int nx,ny,nz;
  ((Patch *)root_patch)->size(&nx,&ny,&nz);
  }
  TRACE3("size = %d %d %d",nx,ny,nz);
}
