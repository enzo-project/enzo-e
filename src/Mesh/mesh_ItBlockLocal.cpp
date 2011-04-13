// $Id: mesh_ItBlockLocal.cpp 1954 2011-01-25 19:54:37Z bordner $
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItBlockLocal.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItBlockLocal
//----------------------------------------------------------------------

#ifndef CONFIG_USE_CHARM

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

ItBlockLocal::ItBlockLocal ( Patch * patch ) throw ()
  : patch_(patch),
    index1_(0)
{}

//----------------------------------------------------------------------

ItBlockLocal::~ItBlockLocal ( ) throw ()
{}

//----------------------------------------------------------------------

Block * ItBlockLocal::operator++ () throw()
{
  index1_ ++;
  if (index1_ > patch_->num_local_blocks()) index1_ = 0;
  return index1_ ? patch_->local_block(index1_ - 1) : 0;
}

//----------------------------------------------------------------------

bool ItBlockLocal::done () const throw()
{
  return index1_ >= patch_->num_local_blocks();
}

#endif
