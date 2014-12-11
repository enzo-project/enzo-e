// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Hierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Implementation of the Hierarchy class

#include "cello.hpp"

#include "mesh.hpp"

#include "charm_mesh.hpp"

//----------------------------------------------------------------------

Hierarchy::Hierarchy 
(
 const Factory * factory,
 int rank, int refinement,
 int process_first, int process_last_plus
 ) throw ()
  :
  factory_((Factory *)factory),
  rank_(rank),
  refinement_(refinement),
  num_blocks_(0),
  block_array_(NULL),
  block_exists_(false),
  block_sync_()
{
  TRACE("Hierarchy::Hierarchy()");
  // Initialize extents
  for (int i=0; i<3; i++) {
    root_size_[i] = 1;
    lower_[i] = 0.0;
    upper_[i] = 1.0;
    blocking_[i] = 0;
  }
}

//----------------------------------------------------------------------

Hierarchy::~Hierarchy() throw()
{
  deallocate_blocks();
}

//----------------------------------------------------------------------

void Hierarchy::pup (PUP::er &p)
{
    
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  const bool up = p.isUnpacking();

  if (up) factory_ = new Factory;
  p | *factory_;
  p | rank_;
  p | refinement_;

  p | num_blocks_;
  // clear if unpacking: load balancing expects num_blocks_ to be
  // updated by CommBlock(CkMigrateMessage) and ~CommBlock(), but
  // checkpoint / restart then double-counts Blocks.

  if (up) num_blocks_ = 0;

  // block_array_ is NULL on non-root processes
  bool allocated=(block_array_ != NULL);
  p|allocated;
  if (allocated) {
    if (up) block_array_=new CProxy_CommBlock;
    p|*block_array_;
  } else {
    block_array_ = NULL;
  }

  p | block_exists_;
  p | block_sync_;
  PUParray(p,root_size_,3);
  PUParray(p,lower_,3);
  PUParray(p,upper_,3);

  PUParray(p,blocking_,3);

}

//----------------------------------------------------------------------

void Hierarchy::set_lower(double x, double y, double z) throw ()
{
  lower_[0] = x;
  lower_[1] = y;
  lower_[2] = z;
}

//----------------------------------------------------------------------

void Hierarchy::set_upper(double x, double y, double z) throw ()
{
  upper_[0] = x;
  upper_[1] = y;
  upper_[2] = z;
}

//----------------------------------------------------------------------

void Hierarchy::set_root_size(int nx, int ny, int nz) throw ()
{
  TRACE3("Hierarchy::set_root_size(%d %d %d)",nx,ny,nz);

  root_size_[0] = nx;
  root_size_[1] = ny;
  root_size_[2] = nz;

}

//----------------------------------------------------------------------

void Hierarchy::set_blocking(int nx, int ny, int nz) throw ()
{
  TRACE3("Hierarchy::set_blocking(%d %d %d)",nx,ny,nz);

  blocking_[0] = nx;
  blocking_[1] = ny;
  blocking_[2] = nz;

}

//----------------------------------------------------------------------

int Hierarchy::rank() const throw ()
{ 
  return rank_; 
}

//----------------------------------------------------------------------

void Hierarchy::root_size(int * nx, int * ny, int * nz) const throw ()
{
  if (nx) *nx = root_size_[0];
  if (ny) *ny = root_size_[1];
  if (nz) *nz = root_size_[2];
}

//----------------------------------------------------------------------

void Hierarchy::lower(double * x, double * y, double * z) const throw ()
{
  if (x) *x = lower_[0];
  if (y) *y = lower_[1];
  if (z) *z = lower_[2];
}
//----------------------------------------------------------------------

void Hierarchy::upper(double * x, double * y, double * z) const throw ()
{
  if (x) *x = upper_[0];
  if (y) *y = upper_[1];
  if (z) *z = upper_[2];
}

//----------------------------------------------------------------------

void Hierarchy::blocking (int * nbx, int * nby, int * nbz) const throw()
{
  if (nbx) (*nbx) = blocking_[0];
  if (nby) (*nby) = blocking_[1];
  if (nbz) (*nbz) = blocking_[2];
}

size_t Hierarchy::num_blocks(int * nbx, 
			     int * nby,
			     int * nbz) const throw()
{ 
  if (nbx) *nbx = blocking_[0];
  if (nby) *nby = blocking_[1];
  if (nbz) *nbz = blocking_[2];

  return blocking_[0]*blocking_[1]*blocking_[2];
}

//----------------------------------------------------------------------

void Hierarchy::deallocate_blocks() throw()
{

  if (block_exists_) {
    block_array_->ckDestroy();
    delete block_array_; block_array_ = 0;
    block_exists_ = false;
  }

}

//----------------------------------------------------------------------

void Hierarchy::create_forest
(
 FieldDescr   * field_descr,
 bool allocate_data,
 bool testing,
 int process_first, int process_last_plus) throw()
{
    // determine block size
    int mbx = root_size_[0] / blocking_[0];
    int mby = root_size_[1] / blocking_[1];
    int mbz = root_size_[2] / blocking_[2];

    // Check that blocks evenly subdivide forest
    if (! ((blocking_[0]*mbx == root_size_[0]) &&
	   (blocking_[1]*mby == root_size_[1]) &&
	   (blocking_[2]*mbz == root_size_[2]))) {

      ERROR6("Forest::allocate_array_()",  
	     "CommBlocks must evenly subdivide forest: "
	     "forest size = (%d %d %d)  block count = (%d %d %d)",
	     root_size_[0],root_size_[1],root_size_[2],
	     blocking_[0],blocking_[1],blocking_[2]);
      
    }

    // CREATE AND INITIALIZE NEW DATA BLOCKS

    int num_field_blocks = 1;

    TRACE("Allocating block_array_");

    block_array_ = new CProxy_CommBlock;

    (*block_array_) = factory_->create_block_array
      (blocking_[0],blocking_[1],blocking_[2],
       mbx,mby,mbz,
       num_field_blocks,
       testing);
    
    block_exists_ = allocate_data;
    block_sync_.set_stop(blocking_[0]*blocking_[1]*blocking_[2]);

}

