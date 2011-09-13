// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Hierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Implementation of the Hierarchy class

#include "cello.hpp"

#include "mesh.hpp"

#include "mesh_charm.hpp"

//----------------------------------------------------------------------

Hierarchy::Hierarchy
(
 const Factory * factory
 ) throw ()
  : factory_(factory),
    tree_(0)
{
  // Initialize extents
  for (int i=0; i<3; i++) {
    lower_[i] = 0.0;
    upper_[i] = 1.0;
  }

}

//----------------------------------------------------------------------

Hierarchy::~Hierarchy() throw()
{
  for (size_t i=0; i<patch_list_.size(); i++) {
    delete patch_list_[i];
    patch_list_[i] = 0;
  }
  delete tree_;
}

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
  return patch_list_.size();
}

//----------------------------------------------------------------------

Patch * Hierarchy::patch(size_t i) throw()
{
  return ( patch_list_.size()-1 >= i )? patch_list_[i] : 0;
}

//----------------------------------------------------------------------

Patch * Hierarchy::patch(size_t i) const throw()
{
  return ( patch_list_.size()-1 >= i )? patch_list_[i] : 0;
}

//----------------------------------------------------------------------

void Hierarchy::create_root_patch 
(
 GroupProcess * group_process,
 FieldDescr   * field_descr,
 int nx, int ny, int nz,
 int nbx, int nby, int nbz
) throw()
{
  //--------------------------------------------------
  // parameters_->group_set (0,"Mesh");
  //--------------------------------------------------
  // parameter: Mesh:refine
  // parameter: Mesh:max_level
  // parameter: Mesh:balanced
  // parameter: Mesh:backfill
  // parameter: Mesh:coalesce
  //--------------------------------------------------

  // mesh_->set_refine_factor (parameters_->value_integer("refine",    2));
  // mesh_->set_max_level     (parameters_->value_integer("max_level", 0));
  // mesh_->set_balanced      (parameters_->value_logical("balanced",  true));
  // mesh_->set_backfill      (parameters_->value_logical("backfill",  true));
  // mesh_->set_coalesce      (parameters_->value_logical("coalesce",  true));

  if (patch_list_.size() != 0) {
    ERROR("Hierarchy::create_root",
	  "Creating new root Patch in non-empty Hierarchy");
  } else {

    Patch * root_patch = factory()->create_patch
      (group_process,
       nx,ny,nz,    // size
       0,0,0,       // offset
       nbx,nby,nbz, // blocking
       lower_[0], lower_[1], lower_[2],
       upper_[0], upper_[1], upper_[2]);
    
    root_patch->allocate_blocks(field_descr);

    insert_patch(root_patch);

    // Parallel layout of the root patch
  
    //--------------------------------------------------
    // parameter: Mesh:root_process_first
    // parameter: Mesh:root_process_count
    //--------------------------------------------------

    // int process_first = parameters_->value_integer("Mesh:root_process_first",0);
    // int process_count = parameters_->value_integer("Mesh:root_process_count",1);

    // patch(0)->layout()->set_process_range(process_first, process_count);

  }
}

//----------------------------------------------------------------------

void Hierarchy::insert_patch(Patch * patch) throw()
{
  int size = patch_list_.size();
  patch_list_.resize(size + 1);
  patch_list_[size] = patch;
}

// //----------------------------------------------------------------------

// bool Hierarchy::balanced() const throw ()
// {
//   return balanced_; 
// }

// //----------------------------------------------------------------------

// void Hierarchy::set_balanced(bool balanced) throw ()
// { 
//   balanced_ = balanced; 
// }

// //----------------------------------------------------------------------

// bool Hierarchy::backfill() const throw ()
// {
//   return backfill_; 
// }

// //----------------------------------------------------------------------

// void Hierarchy::set_backfill(bool backfill) throw ()
// { 
//   backfill_ = backfill; 
// }

// //----------------------------------------------------------------------

// bool Hierarchy::coalesce() const throw ()
// {
//   return coalesce_; 
// }

// //----------------------------------------------------------------------

// void Hierarchy::set_coalesce(bool coalesce) throw ()
// { 
//   coalesce_ = coalesce; 
// }

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
#  include "mesh.def.h"
#endif
