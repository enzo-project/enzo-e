// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItFilePatch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItFilePatch
//----------------------------------------------------------------------

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

ItFilePatch::ItFilePatch ( const Input * input ) throw ()
  : input_(input),
    index1_(0)
{}

//----------------------------------------------------------------------

ItFilePatch::~ItFilePatch ( ) throw ()
{}

//----------------------------------------------------------------------

Patch * ItFilePatch::operator++ () throw()
{
  // index1_ ++;
  // if (index1_ > input_->num_patches()) index1_ = 0;
  // return index1_ ? input_->patch(index1_ - 1) : 0;
}

//----------------------------------------------------------------------

bool ItFilePatch::done () const throw()
{
  // return index1_ >= input_->num_patches();
}




