// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-10
/// @brief    Implementation of the Patch class

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Patch::Patch
(
 const Factory * factory, 
 GroupProcess * group_process,
 int nx,  int ny,  int nz,
 int nx0, int ny0, int nz0,
 int nbx, int nby, int nbz,
 double xm, double ym, double zm,
 double xp, double yp, double zp
) throw()
  :
#ifdef CONFIG_USE_CHARM
  block_exists_(false),
#endif
  factory_(factory),
  group_process_(group_process),
  layout_(new Layout (nbx,nby,nbz))
{

  TRACE0;
 // Check 

  if ( ! ((nx >= nbx) && (ny >= nby) && (nz >= nbz))) {
	     
    ERROR6("Patch::Patch", 
	   "Patch size (%d,%d,%d) must be larger than blocking (%d,%d,%d)",
	   nx,ny,nz,nbx,nby,nbz);
  }

  // create a group process if one doesn't exist
  if (! group_process_) {
    group_process_ = GroupProcess::create();
  }

  // set layout process range

  layout_ -> set_process_range(0,group_process_->size());

  size_[0] = nx;
  size_[1] = ny;
  size_[2] = nz;

  offset_[0] = nx0;
  offset_[1] = ny0;
  offset_[2] = nz0;

  blocking_[0] = nbx;
  blocking_[1] = nby;
  blocking_[2] = nbz;

  lower_[0] = xm;
  lower_[1] = ym;
  lower_[2] = zm;

  upper_[0] = xp;
  upper_[1] = yp;
  upper_[2] = zp;

}

//----------------------------------------------------------------------

Patch::~Patch() throw()
{
  deallocate_blocks();
  delete layout_;
}

//----------------------------------------------------------------------

void Patch::size (int * npx, int * npy, int * npz) const throw()
{
  if (npx) (*npx) = size_[0];
  if (npy) (*npy) = size_[1];
  if (npz) (*npz) = size_[2];
}

//----------------------------------------------------------------------

void Patch::offset (int * nx0, int * ny0, int * nz0) const throw()
{
  if (nx0) (*nx0) = offset_[0];
  if (ny0) (*ny0) = offset_[1];
  if (nz0) (*nz0) = offset_[2];
}

//----------------------------------------------------------------------

void Patch::blocking (int * nbx, int * nby, int * nbz) const throw()
{
  if (nbx) (*nbx) = blocking_[0];
  if (nby) (*nby) = blocking_[1];
  if (nbz) (*nbz) = blocking_[2];
}

//----------------------------------------------------------------------

Layout * Patch::layout () const throw()
{
  return layout_;
}

//----------------------------------------------------------------------
  
void Patch::lower(double * xm, double * ym, double * zm) const throw ()
{
  if (xm) (*xm) = lower_[0];
  if (ym) (*ym) = lower_[1];
  if (zm) (*zm) = lower_[2];
}

//----------------------------------------------------------------------
void Patch::upper(double * xp, double * yp, double * zp) const throw ()
{
  if (xp) (*xp) = upper_[0];
  if (yp) (*yp) = upper_[1];
  if (zp) (*zp) = upper_[2];
}

//----------------------------------------------------------------------
int Patch::index() const throw ()
{
  // ASSUMES ROOT PATCH
  return 0;
}

//======================================================================

void Patch::allocate_array(FieldDescr * field_descr,
			   bool allocate_blocks) throw()
{

#ifndef CONFIG_USE_CHARM

  // determine local block count nb
  int nb = num_local_blocks();

  // create local blocks
  block_.resize(nb);

#endif

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

    ERROR6("Patch::allocate_array",  
	   "Blocks must evenly subdivide Patch: "
	   "patch size = (%d %d %d)  block count = (%d %d %d)",
	   size_[0],size_[1],size_[2],
	   nbx,nby,nbz);
      
  }

  // Determine size of each block
  double xb = (upper_[0] - lower_[0]) / nbx;
  double yb = (upper_[1] - lower_[1]) / nby;
  double zb = (upper_[2] - lower_[2]) / nbz;

  // CREATE AND INITIALIZE NEW DATA BLOCKS

#ifdef CONFIG_USE_CHARM

  // WARNING: assumes patch is on Pe 0!!! (OK for unigrid)

  if (CkMyPe() == 0) {

    block_array_ = factory_->create_block_array
      (nbx,nby,nbz,
       mbx,mby,mbz,
       lower_[0],lower_[1],lower_[2],
       xb,yb,zb,
       allocate_blocks);
    
    block_exists_ = allocate_blocks;

  }

#else

  int ip = group_process_->rank();

  for (int ib=0; ib<nb; ib++) {

    // Get index of block ib in the patch

    int ibx,iby,ibz;
    layout_->block_indices (ip,ib, &ibx, &iby, &ibz);

    // create a new data block

    Block * block = factory_->create_block 
      (ibx,iby,ibz,
       nbx,nby,nbz,
       mbx,mby,mbz,
       lower_[0],lower_[1],lower_[2],
       xb,yb,zb);

    // Store the data block in the block array
    block_[ib] = block;

    // Allocate data on the block

    block->allocate(field_descr);

  }
#endif /* ! CONFIG_USE_CHARM */
}

//----------------------------------------------------------------------

void Patch::deallocate_blocks() throw()
{

#ifdef CONFIG_USE_CHARM

  if (block_exists_) {
    block_array_.ckDestroy();
    block_exists_ = false;
  }

#else

  for (size_t i=0; i<block_.size(); i++) {
    delete block_[i];
    block_[i] = 0;
  }

#endif
}

//----------------------------------------------------------------------
#ifndef CONFIG_USE_CHARM
size_t Patch::num_local_blocks() const  throw()
{
  int rank = group_process_->rank();
  return layout_->local_count(rank);
}
#endif

//----------------------------------------------------------------------

#ifndef CONFIG_USE_CHARM
Block * Patch::local_block(size_t i) const throw()
{
  return (i < block_.size()) ? block_[i] : 0;

}
#endif



