// $Id: mesh_Mesh.cpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Brief description of file mesh_Mesh.cpp

#include "global.hpp"
#include "parameters.hpp"
#include "mesh.hpp"

//----------------------------------------------------------------------

Mesh::Mesh(Global * global,
	   DataDescr * data_descr) throw ()
 : global_(global)
{

  Parameters * parameters = global_->parameters();

  // Initialize the AMR tree

  tree_ = NULL;

  // Initialize dimension_

  parameters->set_current_group ("Physics");

  dimension_ = parameters->value_integer("dimensions");

  // Initialize root_size_[]

  parameters->set_current_group ("Mesh");

  // Patch and block sizes

  min_patch_size_  = parameters->list_value_integer(0,"patch_size",4);
  max_patch_size_  = parameters->list_value_integer(1,"patch_size",128);

  min_block_size_  = parameters->list_value_integer(0,"block_size",4);
  max_block_size_  = parameters->list_value_integer(1,"block_size",128);

  // Mesh size

  root_size_[0] = parameters->list_value_integer(0,"root_size",1);
  root_size_[1] = parameters->list_value_integer(1,"root_size",1);
  root_size_[2] = parameters->list_value_integer(2,"root_size",1);

  // Domain extent for root_patch_

  parameters->set_current_group("Domain");
  double lower[3];
  double upper[3];
  lower[0] = parameters->list_value_scalar(0,"extent",0.0);
  upper[0] = parameters->list_value_scalar(1,"extent",1.0);
  lower[1] = parameters->list_value_scalar(2,"extent",0.0);
  upper[1] = parameters->list_value_scalar(3,"extent",1.0);
  lower[2] = parameters->list_value_scalar(4,"extent",0.0);
  upper[2] = parameters->list_value_scalar(5,"extent",1.0);

  // Create the root patch with the given patchsize and maximum block size

  root_patch_ = new Patch (data_descr,root_size_,max_block_size_,lower,upper);

  // Mesh AMR settings

  refine_    = parameters->value_integer("refine",    2);
  max_level_ = parameters->value_integer("max_level", 0);
  balanced_  = parameters->value_logical("balanced",  true);
  backfill_  = parameters->value_logical("backfill",  true);
  coalesce_  = parameters->value_logical("coalesce",  true);

}

//----------------------------------------------------------------------

Mesh::~Mesh() throw ()
{
}

//----------------------------------------------------------------------

Mesh::Mesh(const Mesh & mesh) throw ()
/// @param     mesh  Object being copied
{
}

//----------------------------------------------------------------------

Mesh & Mesh::operator= (const Mesh & mesh) throw ()
/// @param     mesh  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}
