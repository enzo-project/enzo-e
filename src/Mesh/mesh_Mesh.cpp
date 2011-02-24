// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Brief description of file mesh_Mesh.cpp

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Mesh::Mesh(DataDescr * data_descr) throw ()
  : tree_(0),
    root_patch_(new Patch),
    dimension_(0),
    min_patch_size_(0),
    max_patch_size_(0),
    min_block_size_(0),
    max_block_size_(0),
    refine_(0),
    max_level_(0),
    balanced_(0),
    backfill_(0),
    coalesce_(0)

{
  root_patch_->set_data_descr(data_descr);
  for (int i=0; i<3; i++) {
    lower_[i] = 0.0;
    upper_[i] = 1.0;
    root_size_[i] = 0;
  }
}

//----------------------------------------------------------------------

Mesh::~Mesh() throw()
{
  delete root_patch_;
  delete tree_;
}

