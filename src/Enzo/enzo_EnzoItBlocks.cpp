// $Id: enzo_EnzoItBlocks.cpp 1954 2011-01-25 19:54:37Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoItBlocks.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of EnzoItBlocks

#include "cello.hpp"

#include "enzo.hpp"

EnzoItBlocks::EnzoItBlocks
(
 Patch * patch, 
 EnzoDescr * enzo
 ) throw ()
  : Iterator(), 
    patch_(patch),
    curr_(0),
    enzo_(enzo)
{}

//----------------------------------------------------------------------

EnzoItBlocks::~EnzoItBlocks() throw ()
{
  patch_ = 0; 
  curr_ = 0;
}

//----------------------------------------------------------------------

void * EnzoItBlocks::operator++ ()
{

  if (curr_ == patch_->block_count()) curr_ = 0;

  curr_ ++;

  return patch_->block(curr_ - 1);
}


