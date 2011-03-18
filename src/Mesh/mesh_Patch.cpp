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
 int nbx, int nby, int nbz
) throw()
  : factory_(factory),
    group_process_(group_process),
    layout_(new Layout (nbx,nby,nbz)),
    block_()
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

  lower_[0] = 0.0;
  lower_[1] = 0.0;
  lower_[2] = 0.0;

  upper_[0] = 1.0;
  upper_[1] = 1.0;
  upper_[2] = 1.0;

}

//----------------------------------------------------------------------

Patch::~Patch() throw()
{
  delete layout_;
  deallocate_blocks();
}

//----------------------------------------------------------------------

Patch::Patch(const Patch & patch,
	     FieldDescr *  field_descr) throw()
{
  deallocate_blocks();

  allocate_blocks(field_descr);
}

//----------------------------------------------------------------------

// Patch & Patch::operator= (const Patch & patch) throw()
// {
//   deallocate_blocks();
//   allocate_blocks();
//   return *this;
// }

//----------------------------------------------------------------------

void Patch::size (int * npx, int * npy, int * npz) const throw()
{
  if (npx) *npx = size_[0];
  if (npy) *npy = size_[1];
  if (npz) *npz = size_[2];
}

//----------------------------------------------------------------------

void Patch::blocking (int * nbx, int * nby, int * nbz) const throw()
{
  *nbx = blocking_[0];
  *nby = blocking_[1];
  *nbz = blocking_[2];
}

//----------------------------------------------------------------------

Layout * Patch::layout () const throw()
{
  return layout_;
}

//----------------------------------------------------------------------
  
void Patch::lower(double * xm, double * ym, double * zm) const throw ()
{
  if (xm) *xm = lower_[0];
  if (ym) *ym = lower_[1];
  if (zm) *zm = lower_[2];
}

//----------------------------------------------------------------------
void Patch::set_lower(double xm, double ym, double zm) throw ()
{
  lower_[0] = xm;
  lower_[1] = ym;
  lower_[2] = zm;
};

//----------------------------------------------------------------------
void Patch::upper(double * xp, double * yp, double * zp) const throw ()
{
  if (xp) *xp = upper_[0];
  if (yp) *yp = upper_[1];
  if (zp) *zp = upper_[2];
}

//----------------------------------------------------------------------
void Patch::set_upper(double xp, double yp, double zp) throw ()
{
  upper_[0] = xp;
  upper_[1] = yp;
  upper_[2] = zp;
};

//----------------------------------------------------------------------

void Patch::allocate_blocks(FieldDescr * field_descr) throw()
{

  UNTESTED("Patch::allocate_blocks()");

  // determine local block count nb
  
  int nb = num_blocks();

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
	     "patch size = (%d %d %d)  num_blocks = (%d %d %d)",
	     size_[0],size_[1],size_[2],
	     nbx,nby,nbz);

    ERROR("Patch::allocate_blocks",  buffer);
      
  }

  // Determine size of each block
  double bhx = (upper_[0] - lower_[0]) / nbx;
  double bhy = (upper_[1] - lower_[1]) / nby;
  double bhz = (upper_[2] - lower_[2]) / nbz;

  // CREATE AND INITIALIZE NEW DATA BLOCKS

  for (int ib=0; ib<nb; ib++) {

    // Get index of this block in the patch
    int ibx,iby,ibz;
    layout_->block_indices (ib, &ibx, &iby, &ibz);

    // create a new data block

    Block * block = factory_->create_block (this,field_descr,ibx,iby,ibz,mbx,mby,mbz);

    // Store the data block
    block_[ib] = block;


    // INITIALIZE FIELD BLOCK

    FieldBlock * field_block = block->field_block();

    double xm,xp,ym,yp,zm,zp;

    xm = lower_[0] + ibx*bhx;
    ym = lower_[1] + iby*bhy;
    zm = lower_[2] + ibz*bhz;

    xp = lower_[0] + (ibx+1)*bhx;
    yp = lower_[1] + (iby+1)*bhy;
    zp = lower_[2] + (ibz+1)*bhz;

    block->set_lower(xm,ym,zm);
    block->set_upper(xp,yp,zp);

    // Allocate field data

    field_block->allocate_array();
    field_block->allocate_ghosts();

    WARNING("Patch::allocate_blocks",
		    "allocating all ghosts in patch");

		    
    // INITIALIZE PARTICLE BLOCK

  }
}

//----------------------------------------------------------------------

void Patch::deallocate_blocks() throw()
{
  for (size_t i=0; i<block_.size(); i++) {
    delete block_[i];
    block_[i] = 0;
  }
}

//----------------------------------------------------------------------

bool Patch::blocks_allocated() const throw() 
{
  UNTESTED("Patch::blocks_allocated()");

  bool allocated = true;

  if (block_.size() < num_blocks()) {
      allocated = false;
  } else {
    for (size_t i=0; i<block_.size(); i++) {
      if (block_[i] == NULL) allocated = false;
    }
  }
    
  return allocated;
}

//----------------------------------------------------------------------

size_t Patch::num_blocks() const  throw()
{
  int rank = group_process_->rank();
  return layout_->local_count(rank);
}

//----------------------------------------------------------------------

Block * Patch::block(int i) const throw()
{
  return block_[i];
}

