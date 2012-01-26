// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItBlock
//----------------------------------------------------------------------

#include "mesh.hpp"

//----------------------------------------------------------------------

ItBlock::ItBlock ( Patch * patch ) throw ()
  : patch_(patch)
{}

//----------------------------------------------------------------------

ItBlock::~ItBlock ( ) throw ()
{}

//----------------------------------------------------------------------

Block * ItBlock::operator++ () throw()
{
#ifdef CONFIG_USE_CHARM
  //
  Block * block;
  int nbx,nby,nbz;
  size_t nb = patch_->num_blocks(&nbx,&nby,&nbz);
  do {
    index1_++;
    int ibx = (index1_ - 1) % nbx;
    int iby = (((index1_-1) - ibx)/nbx) % nby;
    int ibz = (index1_-1)/(nbx*nby);
    CProxy_Block block_array = patch_->block_array();
    block = block_array(ibx,iby,ibz).ckLocal();
  } while (block==NULL && index1_ <= nb);
  // assert: (block != NULL) or (index1_ > nb)
  if (index1_ > nb) index1_ = 0;
  // assert: index1_ != 0 implies (block != NULL)
  return index1_ ? block : NULL;
#else
  index1_ ++;
  if (index1_ > patch_->num_local_blocks()) index1_ = 0;
  return index1_ ? patch_->local_block(index1_ - 1) : NULL;
#endif
}

//----------------------------------------------------------------------

bool ItBlock::done () const throw()
{
#ifdef CONFIG_USE_CHARM
  return index1_ >= patch_->num_blocks();
#else
  return index1_ >= patch_->num_local_blocks();
#endif
}

