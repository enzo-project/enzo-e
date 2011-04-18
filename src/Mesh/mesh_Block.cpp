// $Id: mesh_Block.cpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Feb 28 13:22:26 PST 2011
/// @todo     Pre-compute count for p_refresh_face()
/// @todo     Reduce repeated code between p_refresh() and p_refresh_face()
/// @brief    Implementation of the Block object

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Block::Block
(
#ifndef CONFIG_USE_CHARM
 int ix, int iy, int iz,
#endif
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb, // Block width
 int num_field_blocks) throw ()
  : field_block_(),
#ifdef CONFIG_USE_CHARM
    count_refresh_face_(0),
#endif
    cycle_(0),
    time_(0),
    dt_(0)

{ 
  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (nx,ny,nz);
  }

  // Initialize indices into parent patch

#ifdef CONFIG_USE_CHARM
  int ix = thisIndex.x;
  int iy = thisIndex.y;
  int iz = thisIndex.z;
#endif

  index_[0] = ix;
  index_[1] = iy;
  index_[2] = iz;

  // Initialize extent 

  lower_[axis_x] = xpm + ix*xb;
  lower_[axis_y] = ypm + iy*yb;
  lower_[axis_z] = zpm + iz*zb;

  upper_[axis_x] = xpm + (ix+1)*xb;
  upper_[axis_y] = ypm + (iy+1)*yb;
  upper_[axis_z] = zpm + (iz+1)*zb;
}

//----------------------------------------------------------------------

Block::~Block() throw ()
{ 
  // Deallocate field_block_[]
  for (size_t i=0; i<field_block_.size(); i++) {
    delete field_block_[i];
    field_block_[i] = 0;
  }
}

//----------------------------------------------------------------------

Block::Block(const Block & block) throw ()
  : field_block_()
/// @param     block  Object being copied
{
  copy_(block);
}

//----------------------------------------------------------------------

Block & Block::operator = (const Block & block) throw ()
/// @param     block  Source object of the assignment
/// @return    The target assigned object
{
  copy_(block);
  return *this;
}

//----------------------------------------------------------------------

const FieldBlock * Block::field_block (int i) const throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

FieldBlock * Block::field_block (int i) throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

void Block::lower(double * x, double * y, double * z) const throw ()
{
  if (x) *x = lower_[0];
  if (y) *y = lower_[1];
  if (z) *z = lower_[2];
}

//----------------------------------------------------------------------

void Block::upper(double * x, double * y, double * z) const throw ()
{
  if (x) *x = upper_[0];
  if (y) *y = upper_[1];
  if (z) *z = upper_[2];
}

//----------------------------------------------------------------------

// void Block::set_lower(double x, double y, double z) throw ()
// {
//   lower_[0] = x;
//   lower_[1] = y;
//   lower_[2] = z;
// }

// //----------------------------------------------------------------------

// void Block::set_upper(double x, double y, double z) throw ()
// {
//   upper_[0] = x;
//   upper_[1] = y;
//   upper_[2] = z;
// }

//----------------------------------------------------------------------

void Block::index_patch (int * ix=0, int * iy=0, int * iz=0) const throw ()
{
  if (ix) (*ix)=index_[0]; 
  if (iy) (*iy)=index_[1]; 
  if (iz) (*iz)=index_[2]; 
}

//----------------------------------------------------------------------

Block * Block::neighbor (axis_enum axis, face_enum face) const throw()
{
  return NULL;
}

//----------------------------------------------------------------------

void Block::refresh_ghosts(const FieldDescr * field_descr,
			   int index_field_set) throw()
{
  field_block_[index_field_set]->refresh_ghosts(field_descr);
}

//======================================================================
// CHARM FUNCTIONS
//======================================================================

#ifdef CONFIG_USE_CHARM

extern CProxy_Simulation proxy_simulation;
extern CProxy_Main proxy_main;

void Block::p_initial()
{
  TRACE("Block::p_initial");
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  // field_block allocation Was in mesh_Patch() for MPI

  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[0]->allocate_array(field_descr);
    field_block_[0]->allocate_ghosts(field_descr);
  }

  // Apply the initial conditions 

  field_block_[0]->clear(field_descr,1.0);

  Initial * initial = simulation->initial();

  initial->compute(field_descr,this);

  //  field_block_[0]->print(field_descr,"initial");

  initialize(simulation->cycle(), simulation->time());

  // Prepare for the first cycle: perform and disk Output, user
  // Monitoring, apply Stopping criteria [reduction], and compute the
  // next Timestep [reduction]

  // prepare for first cycle: Timestep, Stopping, Monitor, Output

  prepare();
}

//----------------------------------------------------------------------

void Block::prepare()
{

  TRACE("Block::prepare");

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  //--------------------------------------------------
  // Timestep [block]
  //--------------------------------------------------

  double dt_block = simulation->timestep()->compute(field_descr,this);

  //--------------------------------------------------
  // Stopping [block]
  //--------------------------------------------------

  int stop_block = simulation->stopping()->complete(cycle_,time_);

  //--------------------------------------------------
  // Main::p_prepare()
  //--------------------------------------------------

  int num_blocks = simulation->mesh()->patch(0)->num_blocks();

  proxy_main.p_prepare(num_blocks, cycle_, time_, dt_block, stop_block);

}

//----------------------------------------------------------------------

void Block::p_output ()
{
  TRACE("Block::p_output");
}

//----------------------------------------------------------------------

void Block::p_refresh (int nbx, int nby, int nbz, double dt)
{

  TRACE("Block::p_refresh");

  // Update dt_ from Simulation

  dt_ = dt; // (should be updated already?)

  CProxy_Block block_array = thisProxy;

  // Indies of self and neighbors

  int ix = thisIndex.x;
  int iy = thisIndex.y;
  int iz = thisIndex.z;

  int ixm = (ix - 1 + nbx) % nbx;
  int iym = (iy - 1 + nby) % nby;
  int izm = (iz - 1 + nbz) % nbz;
  int ixp = (ix + 1) % nbx;
  int iyp = (iy + 1) % nby;
  int izp = (iz + 1) % nbz;

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  // Get domain lower and upper

  Mesh * mesh = simulation->mesh();

  double lower[3];
  mesh->lower(&lower[0],&lower[1],&lower[2]);
  double upper[3];
  mesh->upper(&upper[0],&upper[1],&upper[2]);

  // Field block size

  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  Boundary * boundary = simulation->boundary();
  const FieldDescr * field_descr = simulation->field_descr();

  FieldFace field_face;

  //--------------------------------------------------
  // X-axis Boundary
  //--------------------------------------------------

  bool periodic = boundary->is_periodic();

  if ( nx > 1) {

    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY

    bool xm_boundary = fabs(lower_[0]-lower[0]) < 1e-7;
    bool xp_boundary = fabs(upper_[0]-upper[0]) < 1e-7;

    // Boundary
    if ( xm_boundary ) boundary->enforce(field_descr,this,face_lower,axis_x);
    if ( xp_boundary ) boundary->enforce(field_descr,this,face_upper,axis_x);

    // Refresh xp <<< xm
    if ( ! xm_boundary || periodic ) {
      field_face.load (field_descr, field_block(), axis_x, face_lower);
      block_array(ixm,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_x, face_upper);
    }

    // Refresh xp >>> xm
    if ( ! xp_boundary || periodic ) {
      field_face.load (field_descr, field_block(), axis_x, face_upper);
      block_array(ixp,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_x, face_lower);
    }

  }

  //--------------------------------------------------
  // Y-axis Boundary
  //--------------------------------------------------

  if ( ny > 1) {

    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY

    bool ym_boundary = fabs(lower_[1]-lower[1]) < 1e-7;
    bool yp_boundary = fabs(upper_[1]-upper[1]) < 1e-7;

    if ( ym_boundary ) boundary->enforce(field_descr, this, face_lower, axis_y);
    if ( yp_boundary ) boundary->enforce(field_descr, this, face_upper, axis_y);

    // Refresh yp <<< ym
    if ( ! ym_boundary || periodic ) {
      field_face.load (field_descr, field_block(), axis_y, face_lower);
      block_array(ix,iym,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_y, face_upper);
    }

    // Refresh yp >>> ym
    if ( ! yp_boundary || periodic ) {
      field_face.load (field_descr, field_block(), axis_y, face_upper);
      block_array(ix,iyp,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_y, face_lower);
    }

  }

  //--------------------------------------------------
  // Z-axis Boundary
  //--------------------------------------------------

  if ( nz > 1) {

    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY

    bool zm_boundary = fabs(lower_[2]-lower[2]) < 1e-7;
    bool zp_boundary = fabs(upper_[2]-upper[2]) < 1e-7;

    if ( zm_boundary ) boundary->enforce(field_descr, this, face_lower, axis_z);
    if ( zp_boundary ) boundary->enforce(field_descr, this, face_upper, axis_z);

    // Refresh zp <<< zm
    if ( ! zm_boundary || periodic ) {
      field_face.load (field_descr, field_block(), axis_z, face_lower);
      block_array(ix,iy,izm).p_refresh_face 
	( field_face.size(), field_face.array(), axis_z, face_upper);
    }

    // Refresh zp >>> zm
    if ( ! zp_boundary || periodic ) {
      field_face.load (field_descr, field_block(), axis_z, face_upper);
      block_array(ix,iy,izp).p_refresh_face 
	(field_face.size(), field_face.array(), axis_z, face_lower);
    }
  }

  // NOTE: p_refresh_face() calls compute, but if no incoming faces
  // it will never get called.  So every block also calls
  // p_refresh_face() itself with a null array

  p_refresh_face (0,0,0,0);

}

//----------------------------------------------------------------------

void Block::p_refresh_face (int n, char * buffer,
			    int axis, int face)
{

  TRACE("Block::p_refresh_face");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  if ( n != 0) {

    // n == 0 is the "default" call from self to ensure p_refresh_face()
    // gets called at least once--required for subsequent compute() call

    FieldFace field_face(n, buffer);

    field_face.store(simulation->field_descr(),
		     field_block(), 
		     axis_enum(axis), 
		     face_enum(face));

  }

  //--------------------------------------------------
  // Count incoming faces
  //--------------------------------------------------

  Mesh * mesh = simulation->mesh();

  double lower[3];
  mesh->lower(&lower[0],&lower[1],&lower[2]);
  double upper[3];
  mesh->upper(&upper[0],&upper[1],&upper[2]);
  
  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  // Determine expected number of incoming Faces
  // (should be function, and possibly precomputed and stored since constant )

  int count = 6 + 1;  // potentially 6 faces + 1 for self

  if (nx==1) count -= 2; // 0D??
  if (ny==1) count -= 2; // 1D
  if (nz==1) count -= 2; // 2D

  bool periodic = simulation->boundary()->is_periodic();

  if (! periodic && nx > 1) {
    if (fabs(lower_[0]-lower[0]) < 1e-7) --count;
    if (fabs(upper_[0]-upper[0]) < 1e-7) --count;
  }
  if (! periodic && ny > 1) {
    if (fabs(lower_[1]-lower[1]) < 1e-7) --count;
    if (fabs(upper_[1]-upper[1]) < 1e-7) --count;
  }
  if (! periodic && nz > 1) {
    if (fabs(lower_[2]-lower[2]) < 1e-7) --count;
    if (fabs(upper_[2]-upper[2]) < 1e-7) --count;
  }

  //--------------------------------------------------
  // Compute
  //--------------------------------------------------

  if (++count_refresh_face_ >= count) {
    count_refresh_face_ = 0;
    compute();
  }
}

//----------------------------------------------------------------------

void Block::compute()
{

  TRACE("Block::compute");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  // DEBUG
  //  FieldDescr * field_descr = simulation->field_descr();
  //  field_block()->print(field_descr,"refresh");

  for (size_t i = 0; i < simulation->num_method(); i++) {

    simulation->method(i) -> compute_block (this,time_,dt_);

  }

  // Update cycle and time

  time_ += dt_;
  ++ cycle_ ;

  // prepare for next cycle: Timestep, Stopping, Monitor, Output

  prepare();

}

//----------------------------------------------------------------------

#endif /* CONFIG_USE_CHARM */

//======================================================================

  void Block::copy_(const Block & block) throw()
{

  // Create a copy of field_block_
  field_block_.resize(block.field_block_.size());
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (*(block.field_block_[i]));
  }
}

//----------------------------------------------------------------------
