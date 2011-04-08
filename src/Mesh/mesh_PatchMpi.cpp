// $Id: mesh_EnzoPatchMpi.cpp 2181 2011-04-07 00:43:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_EnzoPatchMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the EnzoPatchMpi class

#include "cello.hpp"

#include "mesh.hpp"

#ifndef CONFIG_USE_CHARM

//----------------------------------------------------------------------

PatchMpi::PatchMpi
(Factory * factory,
 GroupProcess * group_process,
 int nx,   int ny,  int nz,
 int nbx,  int nby, int nbz,
 double xm, double ym, double zm,
 double xp, double yp, double zp) throw()
  : Patch (factory,group_process,nx,ny,nz,nbx,nby,nbz,xm,ym,zm,xp,yp,zp)
{
}

//----------------------------------------------------------------------

void PatchMpi::allocate_blocks(FieldDescr * field_descr) throw()
{

  // determine local block count nb
  
  int nb = num_local_blocks();

  // create local blocks

  block_.resize(nb);

  // Get number of blocks in the patch
  int nbx,nby,nbz;
  layout_->block_count (&nbx, &nby, &nbz);

  // determine block size
  int mbx = size_[0] / nbx;
  int mby = size_[1] / nby;
  int mbz = size_[2] / nbz;

  // Check that blocks evenly subdivide patch
  if (! ((nbx*mbx == size_[0]) &&
	 (nby*mby == size_[1]) &&
	 (nbz*mbz == size_[2]))) {

    char buffer[ERROR_LENGTH];

    sprintf (buffer,
	     "Blocks must evenly subdivide Patch: "
	     "patch size = (%d %d %d)  block count = (%d %d %d)",
	     size_[0],size_[1],size_[2],
	     nbx,nby,nbz);

    ERROR("PatchMpi::allocate_blocks",  buffer);
      
  }

  // Determine size of each block
  double xb = (upper_[0] - lower_[0]) / nbx;
  double yb = (upper_[1] - lower_[1]) / nby;
  double zb = (upper_[2] - lower_[2]) / nbz;

  // CREATE AND INITIALIZE NEW DATA BLOCKS

  for (int ib=0; ib<nb; ib++) {

    // Get index of this block in the patch
    int ibx,iby,ibz;
    layout_->block_indices (ib, &ibx, &iby, &ibz);

    // create a new data block

    Block * block = factory_->create_block 
      (ibx,iby,ibz,
       mbx,mby,mbz,
       lower_[0],lower_[1],lower_[2],
       xb,yb,zb);

    // Store the data block
    block_[ib] = block;

    // INITIALIZE FIELD BLOCK
    // (move into Block constructor?)
    FieldBlock * field_block = block->field_block();

    // Allocate field data, including ghosts
    
    field_block->allocate_array(field_descr);
    field_block->allocate_ghosts(field_descr);

    // INITIALIZE PARTICLE BLOCK

  }

}

//----------------------------------------------------------------------

PatchMpi::~PatchMpi() throw()
{
  deallocate_blocks();
}

//----------------------------------------------------------------------

void PatchMpi::deallocate_blocks() throw()
{
  for (size_t i=0; i<block_.size(); i++) {
    delete block_[i];
    block_[i] = 0;
  }
}

//----------------------------------------------------------------------

Block * PatchMpi::local_block(size_t i) const throw()
{
  return (i < block_.size()) ? block_[i] : 0;
}

//----------------------------------------------------------------------

#endif /* ! CONFIG_USE_CHARM */

