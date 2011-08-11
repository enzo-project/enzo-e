// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItBlock
//----------------------------------------------------------------------

#ifndef CONFIG_USE_CHARM

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

ItBlock::ItBlock ( Patch * patch ) throw ()
  : patch_(patch),
    index1_(0)
{}

//----------------------------------------------------------------------

ItBlock::~ItBlock ( ) throw ()
{}

//----------------------------------------------------------------------

Block * ItBlock::operator++ () throw()
{
  index1_ ++;
  if (index1_ > patch_->num_local_blocks()) index1_ = 0;
  return index1_ ? patch_->local_block(index1_ - 1) : 0;
}

//----------------------------------------------------------------------

bool ItBlock::done () const throw()
{
  return index1_ >= patch_->num_local_blocks();
}

#endif
