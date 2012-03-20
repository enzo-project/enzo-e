// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Feb 28 13:22:26 PST 2011
/// @brief    Implementation of the Block object
/// @todo     Remove hierarchy dependency via Initial--only need domain bounds

#include "cello.hpp"

#include "mesh.hpp"
#include "main.hpp"

//----------------------------------------------------------------------

Block::Block
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb,    // Block width
 int num_field_blocks) throw ()
  :  num_field_blocks_(num_field_blocks),
     field_block_(),
     cycle_(0),
     time_(0),
#ifdef CONFIG_USE_CHARM
    dt_(0),
    count_refresh_face_(0)
#else
    dt_(0)
#endif

{ 
#ifdef CONFIG_USE_CHARM

#ifdef CONFIG_CHARM_ATSYNC
  usesAtSync = CmiTrue;
#endif

#endif

  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (nx,ny,nz);
  }
  // Initialize indices into parent patch

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

#ifdef CONFIG_USE_CHARM

Block::Block
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb,    // Block width
 int num_field_blocks) throw ()
  : field_block_(),
    num_field_blocks_(num_field_blocks),
    cycle_(0),
    time_(0),
#ifdef CONFIG_USE_CHARM
    dt_(0),
    count_refresh_face_(0)
#else
    dt_(0)
#endif

{ 
#ifdef CONFIG_CHARM_ATSYNC
  usesAtSync = CmiTrue;
#endif

  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (nx,ny,nz);
  }

  // Initialize indices into parent patch

  int ibx = thisIndex.x;
  int iby = thisIndex.y;
  int ibz = thisIndex.z;

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
#endif

//----------------------------------------------------------------------

Block::~Block() throw ()
{ 
  // Deallocate field_block_[]
  for (size_t i=0; i<field_block_.size(); i++) {
    delete field_block_[i];
    field_block_[i] = 0;
  }
  num_field_blocks_ = 0;
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

int Block::index () const throw ()
{
  return index_[0] + size_[0] * (index_[1] + size_[1] * index_[2]);
}

//----------------------------------------------------------------------

void Block::index_patch (int * ix, int * iy, int * iz) const throw ()
{
  if (ix) (*ix) = index_[0]; 
  if (iy) (*iy) = index_[1]; 
  if (iz) (*iz) = index_[2]; 
}

//----------------------------------------------------------------------

void Block::size_patch (int * nx=0, int * ny=0, int * nz=0) const throw ()
{
  if (nx) (*nx)=size_[0]; 
  if (ny) (*ny)=size_[1]; 
  if (nz) (*nz)=size_[2]; 
}

//======================================================================
// MPI FUNCTIONS
//======================================================================

#ifndef CONFIG_USE_CHARM

void Block::refresh_ghosts(const FieldDescr * field_descr,
			   const Patch * patch,
			   int fx, int fy, int fz,
			   int index_field_set) throw()
{
  int ibx,iby,ibz;
  index_patch(&ibx,&iby,&ibz);
  field_block_[index_field_set]
    -> refresh_ghosts (field_descr,
		       patch->group_process(),
		       patch->layout(),
		       ibx,iby,ibz, fx,fy,fz);
}

#endif

//======================================================================
// CHARM FUNCTIONS
//======================================================================

#ifdef CONFIG_USE_CHARM

extern CProxy_Simulation  proxy_simulation;
// extern CProxy_Main        proxy_main;

#endif /* CONFIG_USE_CHARM */

//======================================================================

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_initial()
// dependency: Simulation::field_descr()
// dependency: Simulation: cycle()
// dependency: Simulation: time()
// dependency: Simulation: hierarchy()
//    todo:remove
{
  Simulation * simulation  = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();

  // Initialize the block

  allocate(field_descr);

  // Set the Block cycle and time to match Simulation's

  set_cycle(simulation->cycle());
  set_time (simulation->time());
  set_dt   (simulation->dt());

  // Perform any additional initialization

  initialize ();

  // Apply the initial conditions 

  Initial * initial = simulation->problem()->initial();

  initial->enforce(simulation->hierarchy(),field_descr,this);

  // NOTE: CHARM++ contribute() barrier is to prevent race conditions
  // where Block::p_refresh_face() is called before Block::p_initial()

  // Refresh before prepare()

  contribute( CkCallback(CkIndex_Block::p_call_refresh(), thisProxy) );

}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::prepare()
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  //--------------------------------------------------
  // Enforce boundary conditions
  //--------------------------------------------------

  bool is_boundary[3][2];
  bool axm,axp,aym,ayp,azm,azp;
  determine_boundary_(is_boundary,&axm,&axp,&aym,&ayp,&azm,&azp);
  update_boundary_(is_boundary,axm,axp,aym,ayp,azm,azp);

  FieldDescr * field_descr = simulation->field_descr();

  //--------------------------------------------------
  // Compute local dt
  //--------------------------------------------------

  Problem * problem = simulation->problem();

  double dt_block;
  Timestep * timestep = problem->timestep();

  dt_block = timestep->compute(field_descr,this);

  // Reduce timestep to coincide with scheduled output if needed

  int index_output=0;
  while (Output * output = problem->output(index_output++)) {
    Schedule * schedule = output->schedule();
    dt_block = schedule->update_timestep(time_,dt_block);
  }

  // Reduce timestep to not overshoot final time from stopping criteria

  Stopping * stopping = problem->stopping();

  double time_stop = stopping->stop_time();
  double time_curr = time_;

  dt_block = MIN (dt_block, (time_stop - time_curr));

  //--------------------------------------------------
  // Evaluate local stopping criteria
  //--------------------------------------------------

  int stop_block = stopping->complete(cycle_,time_);

  //--------------------------------------------------
  // Reduce to find Block array minimum dt and stopping criteria
  //--------------------------------------------------

  double min_reduce[2];

  min_reduce[0] = dt_block;
  min_reduce[1] = stop_block ? 1.0 : 0.0;

  CkCallback callback (CkIndex_Block::p_call_output(NULL), thisProxy);
  contribute( 2*sizeof(double), min_reduce, CkReduction::min_double, callback);

}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_call_output(CkReductionMsg * msg)
{

  double * min_reduce = (double * )msg->getData();

  double dt_patch   = min_reduce[0];
  bool   stop_patch = min_reduce[1] == 1.0 ? true : false;

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // WARNING: assumes one patch
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  simulation->update_cycle(cycle_,time_,dt_patch,stop_patch);
 
  // ??? HOW IS cycle_ and time_ update on all processors ensured before index() calls
  // Simulation::p_output()?  Want last block?
  
  // "root" block calls Simulation::p_output()
  if (index() == 0) {
    proxy_simulation.p_output();
  }

  delete msg;
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_refresh (int cycle, double time, double dt)
{
  cycle_ = cycle;
  time_  = time;
  dt_    = dt;

  set_cycle(cycle);
  set_time(time);
  set_dt(dt);

  refresh();
}

//----------------------------------------------------------------------
void Block::p_compute (int cycle, double time, double dt)
{
  cycle_ = cycle;
  time_  = time;
  dt_    = dt;

  set_cycle(cycle);
  set_time(time);
  set_dt(dt);

  compute();
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_call_refresh()
{
  refresh();
}

//--------------------------------------------------

void Block::refresh ()
{

  bool is_boundary[3][2];
  bool axm,axp,aym,ayp,azm,azp;

  determine_boundary_(is_boundary,&axm,&axp,&aym,&ayp,&azm,&azp);

  Simulation * simulation = proxy_simulation.ckLocalBranch();

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

  Boundary * boundary = simulation->problem()->boundary();
  
  bool periodic = boundary->is_periodic();

  CProxy_Block block_array = thisProxy;

  axm = axm && (periodic || ! is_boundary[axis_x][face_lower]);
  axp = axp && (periodic || ! is_boundary[axis_x][face_upper]);
  aym = aym && (periodic || ! is_boundary[axis_y][face_lower]);
  ayp = ayp && (periodic || ! is_boundary[axis_y][face_upper]);
  azm = azm && (periodic || ! is_boundary[axis_z][face_lower]);
  azp = azp && (periodic || ! is_boundary[axis_z][face_upper]);

  // Refresh face ghost zones

  bool lx,ly,lz;
  lx = false;
  ly = false;
  lz = false;

  field_face.set_full(lx,ly,lz);

  FieldDescr * field_descr = simulation->field_descr();

  if (field_descr->refresh_face(2)) {
    if ( axm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, 0, 0);
      block_array(ixm,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), +1, 0, 0);
    }
    if ( axp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, 0, 0);
      block_array(ixp,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), -1, 0, 0);
    }
    if ( aym ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, -1, 0);
      block_array(ix,iym,iz).p_refresh_face 
	(field_face.size(), field_face.array(), 0, +1, 0);
    }
    if ( ayp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, +1, 0);
      block_array(ix,iyp,iz).p_refresh_face 
	(field_face.size(), field_face.array(), 0, -1, 0);
    }
    if ( azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, 0, -1);
      block_array(ix,iy,izm).p_refresh_face 
	(field_face.size(), field_face.array(), 0, 0, +1);
    }
    if ( azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, 0, +1);
      block_array(ix,iy,izp).p_refresh_face 
	(field_face.size(), field_face.array(), 0, 0, -1);
    }
  }

  // Refresh edge ghost zones

  if (field_descr->refresh_face(1)) {
    if ( axm && aym ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, -1, 0);
      block_array(ixm,iym,iz).p_refresh_face 
	(field_face.size(), field_face.array(), +1, +1, 0);
    }
    if ( axm && ayp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, +1, 0);
      block_array(ixm,iyp,iz).p_refresh_face 
	(field_face.size(), field_face.array(), +1, -1, 0);
    }
    if ( axp && aym ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, -1, 0);
      block_array(ixp,iym,iz).p_refresh_face 
	(field_face.size(), field_face.array(), -1, +1, 0);
    }
    if ( axp && ayp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, +1, 0);
      block_array(ixp,iyp,iz).p_refresh_face 
	(field_face.size(), field_face.array(), -1, -1, 0);
    }

    if ( aym && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, -1, -1);
      block_array(ix,iym,izm).p_refresh_face 
	(field_face.size(), field_face.array(), 0, +1, +1);
    }
    if ( aym && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, -1, +1);
      block_array(ix,iym,izp).p_refresh_face 
	(field_face.size(), field_face.array(), 0, +1, -1);
    }
    if ( ayp && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, +1, -1);
      block_array(ix,iyp,izm).p_refresh_face 
	(field_face.size(), field_face.array(), 0, -1, +1);
    }
    if ( ayp && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), 0, +1, +1);
      block_array(ix,iyp,izp).p_refresh_face 
	(field_face.size(), field_face.array(), 0, -1, -1);
    }

    if ( axm && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, 0, -1);
      block_array(ixm,iy,izm).p_refresh_face 
	(field_face.size(), field_face.array(), +1, 0, +1);
    }
    if ( axp && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, 0, -1);
      block_array(ixp,iy,izm).p_refresh_face 
	(field_face.size(), field_face.array(), -1, 0, +1);
    }
    if ( axm && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, 0, +1);
      block_array(ixm,iy,izp).p_refresh_face 
	(field_face.size(), field_face.array(), +1, 0, -1);
    }
    if ( axp && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, 0, +1);
      block_array(ixp,iy,izp).p_refresh_face 
	(field_face.size(), field_face.array(), -1, 0, -1);
    }
  }

  // Refresh corner ghost zones

  if (field_descr->refresh_face(0)) {

    if ( axm && aym && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, -1, -1);
      block_array(ixm,iym,izm).p_refresh_face 
	(field_face.size(), field_face.array(), +1, +1, +1);
    }
    if ( axm && aym && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, -1, +1);
      block_array(ixm,iym,izp).p_refresh_face 
	(field_face.size(), field_face.array(), +1, +1, -1);
    }
    if ( axm && ayp && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, +1, -1);
      block_array(ixm,iyp,izm).p_refresh_face 
	(field_face.size(), field_face.array(), +1, -1, +1);
    }
    if ( axm && ayp && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), -1, +1, +1);
      block_array(ixm,iyp,izp).p_refresh_face 
	(field_face.size(), field_face.array(), +1, -1, -1);
    }
    if ( axp && aym && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, -1, -1);
      block_array(ixp,iym,izm).p_refresh_face 
	(field_face.size(), field_face.array(), -1, +1, +1);
    }
    if ( axp && aym && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, -1, +1);
      block_array(ixp,iym,izp).p_refresh_face 
	(field_face.size(), field_face.array(), -1, +1, -1);
    }
    if ( axp && ayp && azm ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, +1, -1);
      block_array(ixp,iyp,izm).p_refresh_face 
	(field_face.size(), field_face.array(), -1, -1, +1);
    }
    if ( axp && ayp && azp ) {
      FieldFace field_face;
      field_face.load (field_descr, field_block(), +1, +1, +1);
      block_array(ixp,iyp,izp).p_refresh_face 
	(field_face.size(), field_face.array(), -1, -1, -1);
    }
  }

  // NOTE: p_refresh_face() calls compute, but if no incoming faces
  // it will never get called.  So every block also calls
  // p_refresh_face() itself with a null array

  p_refresh_face (0,0,0, 0, 0);

}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::determine_boundary_
(
 bool is_boundary[3][2],
 bool * axm,
 bool * axp,
 bool * aym,
 bool * ayp,
 bool * azm,
 bool * azp
 )
{
  
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  Hierarchy * hierarchy   = simulation->hierarchy();

  double lower_h[3], upper_h[3];
  hierarchy->lower(&lower_h[0],&lower_h[1],&lower_h[2]);
  hierarchy->upper(&upper_h[0],&upper_h[1],&upper_h[2]);

  // return is_boundary[] array of faces on domain boundary

  is_on_boundary(lower_h,upper_h,is_boundary);

  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  // Determine in which directions we need to communicate or update boundary

  *axm = nx > 1;
  *axp = nx > 1;
  *aym = ny > 1;
  *ayp = ny > 1;
  *azm = nz > 1;
  *azp = nz > 1;
}

#endif

//----------------------------------------------------------------------


#ifdef CONFIG_USE_CHARM

void Block::update_boundary_
(
 bool is_boundary[3][2],
 bool axm,
 bool axp,
 bool aym,
 bool ayp,
 bool azm,
 bool azp
)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  Boundary * boundary = simulation->problem()->boundary();
  const FieldDescr * field_descr = simulation->field_descr();


  // Update boundaries

  if ( axm && is_boundary[axis_x][face_lower] ) {
    boundary->enforce(field_descr,this,face_lower,axis_x);
  }
  if ( axp && is_boundary[axis_x][face_upper] ) {
    boundary->enforce(field_descr,this,face_upper,axis_x);
  }
  if ( aym && is_boundary[axis_y][face_lower] ) {
    boundary->enforce(field_descr,this,face_lower,axis_y);
  }
  if ( ayp && is_boundary[axis_y][face_upper] ) {
    boundary->enforce(field_descr,this,face_upper,axis_y);
  }
  if ( azm && is_boundary[axis_z][face_lower] ) {
    boundary->enforce(field_descr,this,face_lower,axis_z);
  }
  if ( azp && is_boundary[axis_z][face_upper] ) {
    boundary->enforce(field_descr,this,face_upper,axis_z);
  }
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_refresh_face (int n, char * buffer, int fx, int fy, int fz)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();

  if ( n != 0) {

    // n == 0 is the call from self to ensure p_refresh_face()
    // always gets called at least once

    bool lx,ly,lz;
    lx = false;
    ly = false;
    lz = false;

    FieldFace field_face(n, buffer);
    field_face.set_full(lx,ly,lz);

    field_face.store (field_descr, field_block(), fx, fy, fz);
  }

  //--------------------------------------------------
  // Count incoming faces
  // (SHOULD NOT RECOMPUTE EVERY CALL)
  //--------------------------------------------------

  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  // Determine axes that may be neighbors

  bool axm = nx > 1;
  bool axp = nx > 1;
  bool aym = ny > 1;
  bool ayp = ny > 1;
  bool azm = nz > 1;
  bool azp = nz > 1;

  // Adjust for boundary faces

  bool periodic = simulation->problem()->boundary()->is_periodic();

  Hierarchy * hierarchy = simulation->hierarchy();

  double lower[3], upper[3];
  hierarchy->lower(&lower[0],&lower[1],&lower[2]);
  hierarchy->upper(&upper[0],&upper[1],&upper[2]);

  bool is_boundary[3][2];
  is_on_boundary (lower,upper,is_boundary);

  axm = axm && (periodic || ! is_boundary[axis_x][face_lower]);
  axp = axp && (periodic || ! is_boundary[axis_x][face_upper]);
  aym = aym && (periodic || ! is_boundary[axis_y][face_lower]);
  ayp = ayp && (periodic || ! is_boundary[axis_y][face_upper]);
  azm = azm && (periodic || ! is_boundary[axis_z][face_lower]);
  azp = azp && (periodic || ! is_boundary[axis_z][face_upper]);

  // Count total expected number of incoming faces

  // self

  int count = 1;

  // faces

  if (field_descr->refresh_face(2)) {
    if ( axm ) ++count;
    if ( axp ) ++count;
    if ( aym ) ++count;
    if ( ayp ) ++count;
    if ( azm ) ++count;
    if ( azp ) ++count;
  }

  // edges

  if (field_descr->refresh_face(1)) {
    if ( axm && aym ) ++count;
    if ( axm && ayp ) ++count;
    if ( axp && aym ) ++count;
    if ( axp && ayp ) ++count;
    if ( aym && azm ) ++count;
    if ( aym && azp ) ++count;
    if ( ayp && azm ) ++count;
    if ( ayp && azp ) ++count;
    if ( azm && axm ) ++count;
    if ( azm && axp ) ++count;
    if ( azp && axm ) ++count;
    if ( azp && axp ) ++count;
  }

  // corners

  if (field_descr->refresh_face(0)) {
    if ( axm && aym && azm ) ++count;
    if ( axm && aym && azp ) ++count;
    if ( axm && ayp && azm ) ++count;
    if ( axm && ayp && azp ) ++count;
    if ( axp && aym && azm ) ++count;
    if ( axp && aym && azp ) ++count;
    if ( axp && ayp && azm ) ++count;
    if ( axp && ayp && azp ) ++count;
  }

  //--------------------------------------------------
  // Compute
  //--------------------------------------------------

  if (++count_refresh_face_ >= count) {
    count_refresh_face_ = 0;
    prepare();
  }
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

// SEE Simulation/simulation_charm_output.cpp for Block::p_write(int)

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::compute()
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

#ifdef CONFIG_USE_PROJECTIONS
  double time_start = CmiWallTimer();
#endif

  FieldDescr * field_descr = simulation->field_descr();

  char buffer[10];
  sprintf (buffer,"%03d-A",cycle_);
  field_block()->print(field_descr,buffer,lower_,upper_);

  int index_method = 0;
  while (Method * method = simulation->problem()->method(index_method++)) {
    method -> compute_block (field_descr,this,time_,dt_);
  }

  sprintf (buffer,"%03d-B",cycle_);
  field_block()->print(field_descr,buffer,lower_,upper_);

#ifdef CONFIG_USE_PROJECTIONS
  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  // Update Block cycle and time

  time_ += dt_;
  ++ cycle_ ;

  set_cycle(cycle_);
  set_time(time_);
  set_dt(dt_);
  
  // prepare for next cycle: Timestep, Stopping, Monitor, Output

  refresh();

}
#endif /* CONFIG_USE_CHARM */

//======================================================================

void Block::copy_(const Block & block) throw()
{

  num_field_blocks_ = block.num_field_blocks_;

  field_block_.resize(block.field_block_.size());
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (*(block.field_block_[i]));
  }

  for (int i=0; i<3; i++) {
    index_[i] = block.index_[i];
    size_[i] = block.size_[i];
    lower_[i] = block.lower_[i];
    upper_[i] = block.upper_[i];
  }

  cycle_ = block.cycle_;
  time_ = block.time_;
  dt_ = block.dt_;

#ifdef CONFIG_USE_CHARM
  count_refresh_face_ = block.count_refresh_face_;
#endif
  
}

//----------------------------------------------------------------------

void Block::is_on_boundary (double lower[3], double upper[3],
			     bool is_boundary[3][2]) throw()
{

  // COMPARISON MAY BE INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY

  is_boundary[axis_x][face_lower] = 
    (cello::err_abs(lower_[axis_x],lower[axis_x]) < 1e-6);
  is_boundary[axis_y][face_lower] = 
    (cello::err_abs(lower_[axis_y],lower[axis_y]) < 1e-6);
  is_boundary[axis_z][face_lower] = 
    (cello::err_abs(lower_[axis_z],lower[axis_z]) < 1e-6);
  is_boundary[axis_x][face_upper] = 
    (cello::err_abs(upper_[axis_x],upper[axis_x]) < 1e-6);
  is_boundary[axis_y][face_upper] = 
    (cello::err_abs(upper_[axis_y],upper[axis_y]) < 1e-6);
  is_boundary[axis_z][face_upper] = 
    (cello::err_abs(upper_[axis_z],upper[axis_z]) < 1e-6);
}
//----------------------------------------------------------------------

void Block::allocate (FieldDescr * field_descr) throw()
{ 
  // Allocate fields
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i]->allocate_array(field_descr);
    field_block_[i]->allocate_ghosts(field_descr);
  }
}

//----------------------------------------------------------------------

