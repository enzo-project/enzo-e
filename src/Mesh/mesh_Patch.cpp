// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @todo     Relocate or remove need for field_block->set_boundary_face() calls
/// @brief    Implementation of the Patch class

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Patch::Patch
(
 DataDescr * data_descr,
 int nx,  int ny,  int nz,
 int nbx, int nby, int nbz
) throw()
  : data_descr_(data_descr),
    layout_(new Layout (nbx,nby,nbz)),
    data_block_()

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

  extents_[0] = 0.0;
  extents_[1] = 1.0;
  extents_[2] = 0.0;
  extents_[3] = 1.0;
  extents_[4] = 0.0;
  extents_[5] = 1.0;

#ifdef CONFIG_USE_MPI
  ip_ = Mpi::rank();
#else
  ip_ = 0;
#endif
}

//----------------------------------------------------------------------

Patch::~Patch() throw()
{
  delete layout_;
  deallocate_blocks();
}

//----------------------------------------------------------------------

Patch::Patch(const Patch & patch) throw()
  :  data_descr_ (patch.data_descr())
{
  deallocate_blocks();

  allocate_blocks();
}

//----------------------------------------------------------------------

Patch & Patch::operator= (const Patch & patch) throw()
{
  deallocate_blocks();
  data_descr_ = patch.data_descr();
  allocate_blocks();
  return *this;
}

//----------------------------------------------------------------------

DataDescr * Patch::data_descr () const throw()
{
  return data_descr_;
}

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

void Patch::set_extents (double xm, double xp,
			 double ym, double yp,
			 double zm, double zp) throw()
{
  extents_[0] = xm;  extents_[1] = xp;
  extents_[2] = ym;  extents_[3] = yp;
  extents_[4] = zm;  extents_[5] = zp;
}

//----------------------------------------------------------------------

void Patch::extents (double * xm, double * xp,
		     double * ym, double * yp,
		     double * zm, double * zp) const throw()
{
  if (xm) *xm = extents_[0];  if (xp) *xp = extents_[1];
  if (ym) *ym = extents_[2];  if (yp) *yp = extents_[3];
  if (zm) *zm = extents_[4];  if (zp) *zp = extents_[5];
}

  
//----------------------------------------------------------------------

void Patch::allocate_blocks() throw()
{

  UNTESTED("Patch::allocate_blocks()");

  // determine local block count nb
  
  int nb = num_blocks();

  // create local blocks

  data_block_.resize(nb);

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
  double bx = (extents_[1] - extents_[0]) / nbx;
  double by = (extents_[3] - extents_[2]) / nby;
  double bz = (extents_[5] - extents_[4]) / nbz;

  // CREATE AND INITIALIZE NEW DATA BLOCKS

  for (int ib=0; ib<nb; ib++) {

    // create a new data block
    DataBlock * data_block = new DataBlock(data_descr_,mbx,mby,mbz);

    // Store the data block
    data_block_[ib] = data_block;

    // Get index of this block in the patch
    int ibx,iby,ibz;
    layout_->block_indices (ib, &ibx, &iby, &ibz);

    // INITIALIZE FIELD BLOCK

    FieldBlock * field_block = data_block->field_block();

    double xm,xp,ym,yp,zm,zp;

    xm = extents_[0] + ibx*bx;
    ym = extents_[2] + iby*by;
    zm = extents_[4] + ibz*bz;

    xp = extents_[0] + (ibx+1)*bx;
    yp = extents_[2] + (iby+1)*by;
    zp = extents_[4] + (ibz+1)*bz;

    field_block->set_extent(xm,xp,ym,yp,zm,zp);


    // Allocate field data

    field_block->allocate_array();
    field_block->allocate_ghosts();

    // Set boundaries

    if (mbx > 1) {
      if (ibx==0)     data_block->set_boundary_face(true,face_lower,axis_x);
      if (ibx==nbx-1) data_block->set_boundary_face(true,face_upper,axis_x);
    }
    if (mby > 1) {
      if (iby==0)     data_block->set_boundary_face(true,face_lower,axis_y);
      if (iby==nby-1) data_block->set_boundary_face(true,face_upper,axis_y);
    }
    if (mbz > 1) {
      if (ibz==0)     data_block->set_boundary_face(true,face_lower,axis_z);
      if (ibz==nbz-1) data_block->set_boundary_face(true,face_upper,axis_z);
    }

    WARNING("Patch::allocate_blocks",
		    "allocating all ghosts in patch");

		    
    // INITIALIZE PARTICLE BLOCK

  }
}

//----------------------------------------------------------------------

void Patch::deallocate_blocks() throw()
{
  for (size_t i=0; i<data_block_.size(); i++) {
    delete data_block_[i];
    data_block_[i] = 0;
  }
}

//----------------------------------------------------------------------

bool Patch::blocks_allocated() const throw() 
{
  UNTESTED("Patch::blocks_allocated()");

  bool allocated = true;

  if (data_block_.size() < num_blocks()) {
      allocated = false;
  } else {
    for (size_t i=0; i<data_block_.size(); i++) {
      if (data_block_[i] == NULL) allocated = false;
    }
  }
    
  return allocated;
}

//----------------------------------------------------------------------

int Patch::num_blocks() const  throw()
{
  return layout_->local_count(ip_);
}

//----------------------------------------------------------------------

DataBlock * Patch::block(int i) const throw()
{
  return data_block_[i];
}

