// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Patch class

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Patch::Patch() throw()
  : data_descr_(0),
    data_block_(),
    layout_(0)

{
  size_[0] = 1;
  size_[1] = 1;
  size_[2] = 1;

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

void Patch::set_data_descr (DataDescr * data_descr) throw()
{
  data_descr_ = data_descr;
}

//----------------------------------------------------------------------

DataDescr * Patch::data_descr () const throw()
{
  return data_descr_;
}

//----------------------------------------------------------------------

void Patch::set_size (int npx, int npy, int npz) throw()
{
  size_[0] = npx;
  size_[1] = npy;
  size_[2] = npz;
}

//----------------------------------------------------------------------

void Patch::size (int * npx, int * npy, int * npz) const throw()
{
  if (npx) *npx = size_[0];
  if (npy) *npy = size_[1];
  if (npz) *npz = size_[2];
}

//----------------------------------------------------------------------

void Patch::set_layout (Layout * layout) throw()
{

  // WARNING: potential for dangling pointer
  layout_ = layout;

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
  extents_[0] = xm;
  extents_[1] = xp;
  extents_[2] = ym;
  extents_[3] = yp;
  extents_[4] = zm;
  extents_[5] = zp;
}

//----------------------------------------------------------------------

void Patch::extents (double * xm, double * xp,
		     double * ym, double * yp,
		     double * zm, double * zp) const throw()
{
  if (xm) *xm = extents_[0];
  if (xp) *xp = extents_[1];
  if (ym) *ym = extents_[2];
  if (yp) *yp = extents_[3];
  if (zm) *zm = extents_[4];
  if (zp) *zp = extents_[5];
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

  FieldDescr * field_descr = data_descr_->field_descr();

  // Determine size of each block
  double bx = (extents_[1] - extents_[0]) / nbx;
  double by = (extents_[3] - extents_[2]) / nby;
  double bz = (extents_[5] - extents_[4]) / nbz;

  // CREATE AND INITIALIZE NEW DATA BLOCKS

  for (int ib=0; ib<nb; ib++) {

    // create a new data block
    DataBlock * data_block = new DataBlock;

    // Store the data block
    data_block_[ib] = data_block;

    // Get index of this block in the patch
    int ibx,iby,ibz;
    layout_->block_indices (ib, &ibx, &iby, &ibz);

    // INITIALIZE FIELD BLOCK

    FieldBlock * field_block = data_block->field_block();

    field_block->set_field_descr(field_descr);

    field_block->set_size(mbx,mby,mbz);

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
      if (ibx==0)     field_block->set_boundary_face(face_lower_x,true);
      if (ibx==nbx-1) field_block->set_boundary_face(face_upper_x,true);
    }
    if (mby > 1) {
      if (iby==0)     field_block->set_boundary_face(face_lower_y,true);
      if (iby==nby-1) field_block->set_boundary_face(face_upper_y,true);
    }
    if (mbz > 1) {
      if (ibz==0)     field_block->set_boundary_face(face_lower_z,true);
      if (ibz==nbz-1) field_block->set_boundary_face(face_upper_z,true);
    }

    WARNING("Patch::allocate_blocks",
		    "allocating all ghosts in patch");

		    
    // INITIALIZE PARTICLE BLOCK

  }
}

//----------------------------------------------------------------------

void Patch::deallocate_blocks() throw()
{
  INCOMPLETE("Patch::deallocate_blocks","");
}

//----------------------------------------------------------------------

bool Patch::blocks_allocated() const throw() 
{
  INCOMPLETE("Patch::blocks_allocated","");
  return false;
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

