// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItPatch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItPatch
//----------------------------------------------------------------------

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

ItPatch::ItPatch ( Hierarchy * hierarchy ) throw ()
  : hierarchy_(hierarchy),
    index1_(0)
{}

//----------------------------------------------------------------------

ItPatch::~ItPatch ( ) throw ()
{}

//----------------------------------------------------------------------

Patch * ItPatch::operator++ () throw()
{
  index1_ ++;
  if (index1_ > hierarchy_->num_patches()) index1_ = 0;
  return index1_ ? hierarchy_->patch(index1_ - 1) : 0;
}

//----------------------------------------------------------------------

bool ItPatch::done () const throw()
{
  return index1_ >= hierarchy_->num_patches();
}




