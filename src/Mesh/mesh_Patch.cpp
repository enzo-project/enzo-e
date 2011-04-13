// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Patch class

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Patch::Patch
(
 Factory * factory, 
 GroupProcess * group_process,
 int nx,  int ny,  int nz,
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
  // Check 
  if ( ! ((nx >= nbx) && (ny >= nby) && (nz >= nbz))) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,
	     "Patch size (%d,%d,%d) must be larger than blocking (%d,%d,%d)",
	     nx,ny,nz,nbx,nby,nbz);
    ERROR("Patch::Patch", buffer);
  }

  size_[0] = nx;
  size_[1] = ny;
  size_[2] = nz;

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

// Patch::Patch(const Patch & patch,
// 	     FieldDescr *  field_descr) throw()
// {
//   deallocate_blocks();

//   allocate_blocks(field_descr);
// }

//----------------------------------------------------------------------

// void Patch::set_lower(double xm, double ym, double zm) throw ()
// {
//   lower_[0] = xm;
//   lower_[1] = ym;
//   lower_[2] = zm;
// };

// //----------------------------------------------------------------------
// void Patch::set_upper(double xp, double yp, double zp) throw ()
// {
//   upper_[0] = xp;
//   upper_[1] = yp;
//   upper_[2] = zp;
// };

// //----------------------------------------------------------------------

// Patch & Patch::operator= (const Patch & patch) throw()
// {
//   deallocate_blocks();
//   allocate_blocks();
//   return *this;
// }

//----------------------------------------------------------------------

void Patch::size (int * npx, int * npy, int * npz) const throw()
{
  if (npx) (*npx) = size_[0];
  if (npy) (*npy) = size_[1];
  if (npz) (*npz) = size_[2];
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

//======================================================================

void Patch::allocate_blocks(FieldDescr * field_descr) throw()
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

    char buffer[ERROR_LENGTH];

    sprintf (buffer,
	     "Blocks must evenly subdivide Patch: "
	     "patch size = (%d %d %d)  block count = (%d %d %d)",
	     size_[0],size_[1],size_[2],
	     nbx,nby,nbz);

    ERROR("Patch::allocate_blocks",  buffer);
      
  }

  // Determine size of each block
  double xb = (upper_[0] - lower_[0]) / nbx;
  double yb = (upper_[1] - lower_[1]) / nby;
  double zb = (upper_[2] - lower_[2]) / nbz;

  // CREATE AND INITIALIZE NEW DATA BLOCKS

#ifdef CONFIG_USE_CHARM

  if (CkMyPe() == 0) {
    TRACE("Allocating block array");
    block_ = CProxy_EnzoBlock::ckNew
      (mbx,mby,mbz,
       lower_[0],lower_[1],lower_[2],
       xb,yb,zb, 1,
       nbx,nby,nbz);
    block_exists_ = true;
  }

#else

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
#endif /* ! CONFIG_USE_CHARM */
}

//----------------------------------------------------------------------

void Patch::deallocate_blocks() throw()
{

#ifdef CONFIG_USE_CHARM

  block_.ckDestroy();
  block_exists_ = false;

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

// #ifdef CONFIG_USE_CHARM
//   if (block_exists_) {
//     return num_blocks();
//   } else {
//     return 0;
//   }
// #else

//----------------------------------------------------------------------

#ifndef CONFIG_USE_CHARM
Block * Patch::local_block(size_t i) const throw()
{
  return (i < block_.size()) ? block_[i] : 0;

}
#endif



