// $Id: mesh_ItBlocks.cpp 1954 2011-01-25 19:54:37Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItBlocks.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItBlocks

#include "cello.hpp"

#include "mesh.hpp"

ItBlocks::ItBlocks ( Patch * patch ) throw ()
  : Iterator(), 
    patch_(patch),
    curr_(0)
{}

//----------------------------------------------------------------------

ItBlocks::~ItBlocks() throw ()
{
  patch_ = 0; 
  curr_ = 0;
}

//----------------------------------------------------------------------

void * ItBlocks::operator++ ()
{
  curr_ ++;

  if (curr_ > patch_->num_blocks()) curr_ = 0;

  return curr_ ? patch_->block(curr_ - 1) : 0;
}


