// $Id: mesh_Block.cpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Feb 28 13:22:26 PST 2011
/// @todo     Pre-compute count for p_refresh_face()
/// @todo     Reduce repeated code between p_refresh() and p_refresh_face()
/// @todo     Remove floating point comparisons for determining boundary faces
/// @todo     Remove need for clearing block values before initial conditions (ghost zones accessed but not initialized)
/// @brief    Implementation of the Block object

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Block::Block
(
#ifndef CONFIG_USE_CHARM
 int ibx, int iby, int ibz,
#endif
 int nbx, int nby, int nbz,
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
#ifndef CONFIG_LOAD_BALANCE
  usesAtSync = CmiTrue;
#endif

  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (nx,ny,nz);
  }

  // Initialize indices into parent patch

#ifdef CONFIG_USE_CHARM
  int ibx = thisIndex.x;
  int iby = thisIndex.y;
  int ibz = thisIndex.z;
#endif

  size_[0] = nbx;
  size_[1] = nby;
  size_[2] = nbz;

  index_[0] = ibx;
  index_[1] = iby;
  index_[2] = ibz;

  // Initialize extent 

  lower_[axis_x] = xpm + ibx*xb;
  lower_[axis_y] = ypm + iby*yb;
  lower_[axis_z] = zpm + ibz*zb;

  upper_[axis_x] = xpm + (ibx+1)*xb;
  upper_[axis_y] = ypm + (iby+1)*yb;
  upper_[axis_z] = zpm + (ibz+1)*zb;
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

void Block::size_patch (int * nx=0, int * ny=0, int * nz=0) const throw ()
{
  if (nx) (*nx)=size_[0]; 
  if (ny) (*ny)=size_[1]; 
  if (nz) (*nz)=size_[2]; 
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

  // SHOULD NOT NEED THIS
  field_block_[0]->clear(field_descr,0.001);

  simulation->initial()->compute(field_descr,this);

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

  // Reduce timestep to coincide with scheduled output if needed

  for (size_t i=0; i<simulation->num_output(); i++) {
    Output * output = simulation->output(i);
    dt_block = output->update_timestep(time_,dt_block);
  }

  // Reduce timestep to coincide with end of simulation if needed

  dt_block = MIN (dt_block, (simulation->stopping()->stop_time() - time_));

  //--------------------------------------------------
  // Stopping [block]
  //--------------------------------------------------

  int stop_block = simulation->stopping()->complete(cycle_,time_);
  // DEBUG
  if (stop_block) {
    FieldDescr * field_descr = simulation->field_descr();
  //   field_block()->print(field_descr,"compute",lower_,upper_);
    field_block()->image(field_descr,"final",cycle_,
			 thisIndex.x,thisIndex.y,thisIndex.z);
  }


  //--------------------------------------------------
  // Main::p_prepare()
  //--------------------------------------------------

  int num_blocks = simulation->mesh()->patch(0)->num_blocks();


#undef CHARM_REDUCTION

#ifdef CHARM_REDUCTION
  printf ("contribute(%g)\n",dt_block);
  contribute (sizeof(double),&dt_block,CkReduction::min_double);
#else
  proxy_main.p_prepare(num_blocks, cycle_, time_, dt_block, stop_block);
#endif


}

//----------------------------------------------------------------------

void Block::p_refresh (double dt)
{

  TRACE("Block::p_refresh");

  // Update dt_ from Simulation

  dt_ = dt; // (should be updated already?)

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  //--------------------------------------------------
  // Boundary
  //--------------------------------------------------

  bool boundary_face[3][2];
  
  boundary_face[axis_x][face_lower] = false;
  boundary_face[axis_x][face_upper] = false;
  boundary_face[axis_y][face_lower] = false;
  boundary_face[axis_y][face_upper] = false;
  boundary_face[axis_z][face_lower] = false;
  boundary_face[axis_z][face_upper] = false;

  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  Mesh *             mesh        = simulation->mesh();
  Boundary *         boundary    = simulation->boundary();
  const FieldDescr * field_descr = simulation->field_descr();

  double lower[3];
  mesh->lower(&lower[0],&lower[1],&lower[2]);
  double upper[3];
  mesh->upper(&upper[0],&upper[1],&upper[2]);


  if ( nx > 1) {
    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY
    boundary_face[axis_x][face_lower] = fabs(lower_[0]-lower[0]) < 1e-7;
    boundary_face[axis_x][face_upper] = fabs(upper_[0]-upper[0]) < 1e-7;
    if ( boundary_face[axis_x][face_lower] ) boundary->enforce(field_descr,this,face_lower,axis_x);
    if ( boundary_face[axis_x][face_upper] ) boundary->enforce(field_descr,this,face_upper,axis_x);
  }
  if ( ny > 1) {
    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY
    boundary_face[axis_y][face_lower] = fabs(lower_[1]-lower[1]) < 1e-7;
    boundary_face[axis_y][face_upper] = fabs(upper_[1]-upper[1]) < 1e-7;
    if ( boundary_face[axis_y][face_lower] ) boundary->enforce(field_descr,this,face_lower,axis_y);
    if ( boundary_face[axis_y][face_upper] ) boundary->enforce(field_descr,this,face_upper,axis_y);
  }
  if ( nz > 1) {
    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY
    boundary_face[axis_z][face_lower] = fabs(lower_[2]-lower[2]) < 1e-7;
    boundary_face[axis_z][face_upper] = fabs(upper_[2]-upper[2]) < 1e-7;
    if ( boundary_face[axis_z][face_lower] ) boundary->enforce(field_descr,this,face_lower,axis_z);
    if ( boundary_face[axis_z][face_upper] ) boundary->enforce(field_descr,this,face_upper,axis_z);
  }

  //--------------------------------------------------
  // Refresh
  //--------------------------------------------------

  int ix = thisIndex.x;
  int iy = thisIndex.y;
  int iz = thisIndex.z;

  int nbx = size_[0];
  int nby = size_[1];
  int nbz = size_[2];
  
  int ixm = (ix - 1 + nbx) % nbx;
  int iym = (iy - 1 + nby) % nby;
  int izm = (iz - 1 + nbz) % nbz;
  int ixp = (ix + 1) % nbx;
  int iyp = (iy + 1) % nby;
  int izp = (iz + 1) % nbz;

  FieldFace field_face;

  bool periodic = boundary->is_periodic();

  CProxy_Block block_array = thisProxy;

  if ( nx > 1) {
    // xp <<< xm
    if ( ! boundary_face[axis_x][face_lower] || periodic ) {
      field_face.load (field_descr, field_block(), axis_x, face_lower);
      block_array(ixm,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_x, face_upper);
    }
    // xp >>> xm
    if ( ! boundary_face[axis_x][face_upper] || periodic ) {
      field_face.load (field_descr, field_block(), axis_x, face_upper);
      block_array(ixp,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_x, face_lower);
    }
  }
  if ( ny > 1) {
    // yp <<< ym
    if ( ! boundary_face[axis_y][face_lower] || periodic ) {
      field_face.load (field_descr, field_block(), axis_y, face_lower);
      block_array(ix,iym,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_y, face_upper);
    }
    // yp >>> ym
    if ( ! boundary_face[axis_y][face_upper] || periodic ) {
      field_face.load (field_descr, field_block(), axis_y, face_upper);
      block_array(ix,iyp,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_y, face_lower);
    }
  }
  if ( nz > 1) {
    // zp <<< zm
    if ( ! boundary_face[axis_z][face_lower] || periodic ) {
      field_face.load (field_descr, field_block(), axis_z, face_lower);
      block_array(ix,iy,izm).p_refresh_face 
	( field_face.size(), field_face.array(), axis_z, face_upper);
    }
    // zp >>> zm
    if ( ! boundary_face[axis_z][face_upper] || periodic ) {
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

void Block::p_output_accum (int index_output)
{
  char buffer[80];
  sprintf (buffer,"Block::p_output_accum(%d)",index_output);
  INCOMPLETE(buffer);

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  simulation->output(index_output)->accum_block(this);

  // Synchronize via main chare before writing
  int num_blocks = simulation->mesh()->patch(0)->num_blocks();
  proxy_main.p_output_reduce (num_blocks);
}

//----------------------------------------------------------------------

void Block::compute()
{

  TRACE("Block::compute");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

#ifdef CONFIG_USE_PROJECTIONS
  double time_start = CmiWallTimer();
#endif

  for (size_t i = 0; i < simulation->num_method(); i++) {

    simulation->method(i) -> compute_block (this,time_,dt_);

  }

#ifdef CONFIG_USE_PROJECTIONS
  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif
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
