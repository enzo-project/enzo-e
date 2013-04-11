// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 18:06:42 PST 2011
/// @brief    Implementation of ItBlock
//----------------------------------------------------------------------

#include "mesh.hpp"

//----------------------------------------------------------------------

#ifdef REMOVE_PATCH
ItBlock::ItBlock ( const Hierarchy * hierarchy ) throw ()
  : hierarchy_((Hierarchy * )hierarchy)
#else /* REMOVE_PATCH */
ItBlock::ItBlock ( const Patch * patch ) throw ()
  : patch_((Patch * )patch)
#endif /* REMOVE_PATCH */

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
#ifdef REMOVE_PATCH
  hierarchy_->root_size(&nbx,&nby,&nbz);
  size_t nb = nbx*nby*nbz;
#else /* REMOVE_PATCH */
  size_t nb = patch_->num_blocks(&nbx,&nby,&nbz);
#endif /* REMOVE_PATCH */

  do {
    index1_++;
    int ibx = (index1_ - 1) % nbx;
    int iby = (((index1_-1) - ibx)/nbx) % nby;
    int ibz = (index1_-1)/(nbx*nby);
#ifdef REMOVE_PATCH
    CProxy_CommBlock * block_array = hierarchy_->block_array();
#else /* REMOVE_PATCH */
    CProxy_CommBlock * block_array = patch_->block_array();
#endif /* REMOVE_PATCH */

    block = (*block_array)(ibx,iby,ibz).ckLocal();
  } while (block==NULL && index1_ <= nb);
  // assert: (block != NULL) or (index1_ > nb)
  if (index1_ > nb) index1_ = 0;
  // assert: index1_ != 0 implies (block != NULL)
  return index1_ ? block : NULL;

#else /* CONFIG_USE_CHARM */

#ifdef REMOVE_PATCH
  index1_ ++;
  int nb = hierarchy_->num_local_blocks();
  if (index1_ > nb) index1_ = 0;
  return index1_ ? hierarchy_->local_block(index1_ - 1) : NULL;
#else /* REMOVE_PATCH */
  index1_ ++;
  if (index1_ > patch_->num_local_blocks()) index1_ = 0;
  return index1_ ? patch_->local_block(index1_ - 1) : NULL;
#endif

#endif /* CONFIG_USE_CHARM */
}

//----------------------------------------------------------------------

bool ItBlock::done () const throw()
{
#ifdef REMOVE_PATCH

  int nbx,nby,nbz;
  hierarchy_->root_size(&nbx,&nby,&nbz);
  size_t nb = nbx*nby*nbz;
  return index1_ >= nb;

#else /* REMOVE_PATCH */

# ifdef CONFIG_USE_CHARM
  return index1_ >= patch_->num_blocks();
# else
  return index1_ >= patch_->num_local_blocks();
# endif

#endif /* REMOVE_PATCH */
}

