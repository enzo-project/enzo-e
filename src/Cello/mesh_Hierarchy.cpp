// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Hierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Implementation of the Hierarchy class

#include "cello.hpp"

#include "mesh.hpp"

#include "charm_mesh.hpp"

// #define CELLO_TRACE

//----------------------------------------------------------------------

Hierarchy::Hierarchy 
(
 const Factory * factory,
 int rank, int refinement,
 int max_level) throw ()
  :
  factory_((Factory *)factory),
  rank_(rank),
  refinement_(refinement),
  max_level_(max_level),
  num_blocks_(0),
  num_blocks_level_(),
  num_particles_(0),
  num_zones_total_(0),
  num_zones_real_(0),
  block_array_(),
  block_exists_(false)
{
  TRACE("Hierarchy::Hierarchy()");
  // Initialize extents
				   
  num_blocks_level_.resize(max_level+1);

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
  p | max_level_;

  p | num_blocks_;
  p | num_blocks_level_;
  p | num_particles_;
  p | num_zones_total_;
  p | num_zones_real_;

  // clear if unpacking: load balancing expects num_blocks_ to be
  // updated by Block(CkMigrateMessage) and ~Block(), but
  // checkpoint / restart then double-counts Blocks.

  if (up) {
    for (int i=0; i<num_blocks_level_.size(); i++)
      num_blocks_level_[i]=0;
    num_blocks_      = 0;
    num_particles_   = 0;
    num_zones_total_ = 0;
    num_zones_real_  = 0;
  }

  p | block_array_;
  p | block_exists_;

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

void Hierarchy::root_blocks (int * nbx, int * nby, int * nbz) const throw()
{
  if (nbx) (*nbx) = blocking_[0];
  if (nby) (*nby) = blocking_[1];
  if (nbz) (*nbz) = blocking_[2];
}

//----------------------------------------------------------------------

void Hierarchy::deallocate_blocks() throw()
{
}

//----------------------------------------------------------------------

#ifdef NEW_MSG_REFINE
CProxy_Block Hierarchy::new_block_proxy ( bool allocate_data) throw()
{
  TRACE("Creating block_array_");

  DataMsg * data_msg = NULL;

  block_array_ = factory_->new_block_proxy
    ( data_msg,  blocking_[0],blocking_[1],blocking_[2]);
  return block_array_;
}
#endif

//----------------------------------------------------------------------

void Hierarchy::create_block_array ( bool allocate_data) throw()
{
  // determine block size
  const int mbx = root_size_[0] / blocking_[0];
  const int mby = root_size_[1] / blocking_[1];
  const int mbz = root_size_[2] / blocking_[2];

  // Check that blocks evenly subdivide array
  if (! ((blocking_[0]*mbx == root_size_[0]) &&
	 (blocking_[1]*mby == root_size_[1]) &&
	 (blocking_[2]*mbz == root_size_[2]))) {

    ERROR6("Hierarchy::create_block_array()",  
	   "Blocks must evenly subdivide array: "
	   "array size = (%d %d %d)  block count = (%d %d %d)",
	   root_size_[0],
	   root_size_[1],
	   root_size_[2],
	   blocking_[0],
	   blocking_[1],
	   blocking_[2]);
  }

  // CREATE AND INITIALIZE NEW DATA BLOCKS

  int num_field_blocks = 1;

  TRACE("Allocating block_array_");

  DataMsg * data_msg = NULL;

#ifdef NEW_MSG_REFINE  
  factory_->create_block_array
    ( data_msg,
      block_array_,
      blocking_[0],blocking_[1],blocking_[2],
      mbx,mby,mbz,
      num_field_blocks);
#else
  block_array_ = factory_->create_block_array
    ( data_msg,
      blocking_[0],blocking_[1],blocking_[2],
      mbx,mby,mbz,
      num_field_blocks);
#endif

  block_exists_ = allocate_data;
}

//----------------------------------------------------------------------

void Hierarchy::create_subblock_array
(bool allocate_data, int min_level) throw()
{
  // determine block size

  const int mbx = root_size_[0] / blocking_[0];
  const int mby = root_size_[1] / blocking_[1];
  const int mbz = root_size_[2] / blocking_[2];

  // CREATE AND INITIALIZE NEW DATA BLOCKS

  int num_field_blocks = 1;

  TRACE("Allocating sub-block_array_");

  DataMsg * data_msg = NULL;

  factory_->create_subblock_array
    (data_msg,
     block_array_,min_level,
     blocking_[0],blocking_[1],blocking_[2],
     mbx,mby,mbz,
     num_field_blocks);
    
}

