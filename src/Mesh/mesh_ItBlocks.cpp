// $Id: mesh_ItBlocks.cpp 1954 2011-01-25 19:54:37Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItBlocks.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItBlocks

#include "cello.hpp"

#include "mesh.hpp"

ItBlocks::ItBlocks ( Patch * patch ) throw ()
  : patch_(patch),
    index1_(0)
{}

//----------------------------------------------------------------------

ItBlocks::~ItBlocks() throw ()
{
  patch_ = 0; 
  index1_ = 0;
}

//----------------------------------------------------------------------

void ItBlocks::first () throw()
{
  index1_ = 0;
}

//----------------------------------------------------------------------

void ItBlocks::next () throw()
{
  ++index1_;
}

//----------------------------------------------------------------------

bool ItBlocks::done () const throw()
{
  return index1_ >= patch_->num_blocks();
}

//----------------------------------------------------------------------

DataBlock * ItBlocks::curr () throw()
{
  return (index1_ < patch_->num_blocks()) ? patch_->block(index1_) : 0;
}

//----------------------------------------------------------------------

const DataBlock * ItBlocks::curr () const throw()
{
  return (index1_ < patch_->num_blocks()) ? patch_->block(index1_) : 0;
}

//----------------------------------------------------------------------
int ItBlocks::index (int * ibx, int * iby, int * ibz) throw()
{
  if (index1_) {
      // need to implement mapping of local block 1D index in patch 
      // to global 3D index.  Either use Patch or ItBlock::index()
    INCOMPLETE("ItBlocks::index","");
    //    patch_->layout()->block_indices
  } else {
    WARNING ("ItBlocks::index","Trying to get index of the 'null DataBlock'");
  }
  return 0;
}



