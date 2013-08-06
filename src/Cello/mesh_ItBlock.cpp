// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItBlock
//----------------------------------------------------------------------

#include "mesh.hpp"

//----------------------------------------------------------------------

ItBlock::ItBlock ( const Hierarchy * hierarchy ) throw ()
  : hierarchy_((Hierarchy * )hierarchy)
{}

//----------------------------------------------------------------------

ItBlock::~ItBlock ( ) throw ()
{}

//----------------------------------------------------------------------

CommBlock * ItBlock::operator++ () throw()
{
#ifdef CONFIG_USE_CHARM
  //
  CommBlock * block;
  int nbx,nby,nbz;
  hierarchy_->root_size(&nbx,&nby,&nbz);
  size_t nb = nbx*nby*nbz;

  do {
    index1_++;
    int ibx = (index1_ - 1) % nbx;
    int iby = (((index1_-1) - ibx)/nbx) % nby;
    int ibz = (index1_-1)/(nbx*nby);
    CProxy_CommBlock * block_array = hierarchy_->block_array();

    Index index(ibx,iby,ibz);
    block = (*block_array)[index].ckLocal();

  } while (block==NULL && index1_ <= nb);
  // assert: (block != NULL) or (index1_ > nb)
  if (index1_ > nb) index1_ = 0;
  // assert: index1_ != 0 implies (block != NULL)
  return index1_ ? block : NULL;

#else /* CONFIG_USE_CHARM */

  index1_ ++;
  size_t nb = hierarchy_->num_local_blocks();
  if (index1_ > nb) index1_ = 0;
  return index1_ ? hierarchy_->local_block(index1_ - 1) : NULL;

#endif /* CONFIG_USE_CHARM */
}

//----------------------------------------------------------------------

bool ItBlock::done () const throw()
{
  int nbx,nby,nbz;
  hierarchy_->root_size(&nbx,&nby,&nbz);
  size_t nb = nbx*nby*nbz;
  return index1_ >= nb;
}

