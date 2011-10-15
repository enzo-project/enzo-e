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
#include "main.hpp"

//----------------------------------------------------------------------

#define TEMP_CLEAR_VALUE std::numeric_limits<float>::min() /* in field_FieldBlock.cpp and  mesh_Block.cpp */

//----------------------------------------------------------------------

Block::Block
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb,    // Block width
 int num_field_blocks) throw ()
  :  field_block_(),
#ifdef CONFIG_USE_CHARM
     num_field_blocks_(num_field_blocks),
     count_refresh_face_(0),
     count_refresh_face_x_(0),
     count_refresh_face_y_(0),
     count_refresh_face_z_(0),
#endif
     cycle_(0),
     time_(0),
     dt_(0)

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

#ifdef CONFIG_USE_CHARM
  // WARNING: this constructor should only be called by test code
  ibx = 0;
  iby = 0;
  ibz = 0;
#endif

  size_[0] = nbx;
  size_[1] = nby;
  size_[2] = nbz;

#ifndef CONFIG_USE_CHARM
  index_[0] = ibx;
  index_[1] = iby;
  index_[2] = ibz;
#endif

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
    count_refresh_face_(0),
    count_refresh_face_x_(0),
    count_refresh_face_y_(0),
    count_refresh_face_z_(0),
    cycle_(0),
    time_(0),
    dt_(0)

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

#ifndef CONFIG_USE_CHARM
  index_[0] = ibx;
  index_[1] = iby;
  index_[2] = ibz;
#endif

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
#ifdef CONFIG_USE_CHARM
  num_field_blocks_ = 0;
#endif
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

int Block::index () const throw ()
{
#ifdef CONFIG_USE_CHARM
  return thisIndex.x + size_[0] * (thisIndex.y + size_[1] * thisIndex.z);
#else
  return index_[0] + size_[0] * (index_[1] + size_[1] * index_[2]);
#endif
}

//----------------------------------------------------------------------

void Block::index_patch (int * ix, int * iy, int * iz) const throw ()
{
#ifdef CONFIG_USE_CHARM
  if (ix) (*ix) = thisIndex.x;
  if (iy) (*iy) = thisIndex.y;
  if (iz) (*iz) = thisIndex.z;
#else
  if (ix) (*ix) = index_[0]; 
  if (iy) (*iy) = index_[1]; 
  if (iz) (*iz) = index_[2]; 
#endif
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

//======================================================================
// MPI FUNCTIONS
//======================================================================

#ifndef CONFIG_USE_CHARM

void Block::refresh_ghosts(const FieldDescr * field_descr,
			   const Patch * patch,
			   face_enum face,
			   axis_enum axis,
			   int index_field_set) throw()
{
  int ibx,iby,ibz;
  index_patch(&ibx,&iby,&ibz);
  field_block_[index_field_set]
    -> refresh_ghosts (field_descr,patch, ibx,iby,ibz, face,axis);
}

#endif

//======================================================================
// CHARM FUNCTIONS
//======================================================================

#ifdef CONFIG_USE_CHARM

extern CProxy_Simulation  proxy_simulation;
extern CProxy_Main        proxy_main;

#endif /* CONFIG_USE_CHARM */

//======================================================================

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_initial()
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  // Allocate FieldBlock arrays

  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i]->allocate_array(field_descr);
    field_block_[i]->allocate_ghosts(field_descr);
    WARNING1("Block::p_initial","Clearing field block values to %g",TEMP_CLEAR_VALUE);
    field_block_[i]->clear(field_descr,TEMP_CLEAR_VALUE);
  }

  // Initialize the block

  initialize(simulation->cycle(), simulation->time());

  // Apply the initial conditions 

  Initial * initial = simulation->initial();

  initial->compute(field_descr,this);

#ifdef ORIGINAL_REFRESH  

  // Prepare for first cycle: Timestep, Stopping, Monitor, Output

  prepare();

#else

  // Refresh before prepare()
  axis_enum axis_set = (simulation->temp_update_all()) ? axis_all : axis_x;
  refresh(axis_set);

#endif
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::prepare()
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  //--------------------------------------------------
  // Timestep [block]
  //--------------------------------------------------

  double dt_block;
  dt_block = simulation->timestep()->compute(field_descr,this);

  // Reduce timestep to coincide with scheduled output if needed

  for (size_t i=0; i<simulation->num_output(); i++) {
    Schedule * schedule = simulation->output(i)->schedule();
    dt_block = schedule->update_timestep(time_,dt_block);
  }

  // Reduce timestep to not overshoot final time from stopping criteria

  double time_stop = simulation->stopping()->stop_time();
  double time_curr = time_;

  dt_block = MIN (dt_block, (time_stop - time_curr));

  //--------------------------------------------------
  // Stopping [block]
  //--------------------------------------------------

  int stop_block = simulation->stopping()->complete(cycle_,time_);

  //--------------------------------------------------
  // Main::p_prepare()
  //--------------------------------------------------

  int num_blocks = simulation->hierarchy()->patch(0)->num_blocks();

  simulation->proxy_block_reduce().p_prepare
    (num_blocks, cycle_, time_, dt_block, stop_block);

}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::refresh_axis (axis_enum axis)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  //--------------------------------------------------
  // Boundary
  //--------------------------------------------------

  // WARNING similar code to EnzoSimulationMpi::run() for MPI

  bool boundary_face[2];
  
  boundary_face[face_lower] = false;
  boundary_face[face_upper] = false;

  int n3[3];
  field_block()->size (&n3[0],&n3[1],&n3[2]);

  Hierarchy *        hierarchy   = simulation->hierarchy();
  Boundary *         boundary    = simulation->boundary();
  const FieldDescr * field_descr = simulation->field_descr();

  double lower[3];
  hierarchy->lower(&lower[0],&lower[1],&lower[2]);
  double upper[3];
  hierarchy->upper(&upper[0],&upper[1],&upper[2]);

  bool is_active = n3[axis] > 1;

  if ( is_active ) {
    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY
    boundary_face[face_lower] = fabs(lower_[axis]-lower[axis]) < 1e-7;
    boundary_face[face_upper] = fabs(upper_[axis]-upper[axis]) < 1e-7;
    if ( boundary_face[face_lower] ) boundary->enforce(field_descr,this,face_lower,axis);
    if ( boundary_face[face_upper] ) boundary->enforce(field_descr,this,face_upper,axis);
  }

  //--------------------------------------------------
  // Refresh
  //--------------------------------------------------

  int i3[3]  = { thisIndex.x, thisIndex.y, thisIndex.z};
  int im3[3] = { thisIndex.x, thisIndex.y, thisIndex.z};
  int ip3[3] = { thisIndex.x, thisIndex.y, thisIndex.z};
  int nb3[3] = { size_[0], size_[1], size_[2]};

  im3[axis] = (i3[axis] - 1 + nb3[axis]) % nb3[axis];
  ip3[axis] = (i3[axis] + 1 + nb3[axis]) % nb3[axis];

  FieldFace field_face;

  bool periodic = boundary->is_periodic();

  CProxy_Block block_array = thisProxy;

  bool update_full = simulation->temp_update_full();

  if ( is_active ) {
    if ( ! boundary_face[face_lower] || periodic ) {
      field_face.load (field_descr, field_block(), axis, face_lower, 
		       update_full,update_full);
      block_array(im3[0],im3[1],im3[2]).p_refresh_face 
	(field_face.size(), field_face.array(), axis, face_upper, axis);

    }
    if ( ! boundary_face[face_upper] || periodic ) {
      field_face.load (field_descr, field_block(), axis, face_upper, 
		       update_full,update_full);
      block_array(ip3[0],ip3[1],ip3[2]).p_refresh_face 
	(field_face.size(), field_face.array(), axis, face_lower, axis);
    }
  }
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_compute (double dt, int axis_set)
{
  if (dt != -1) dt_ = dt;
  compute(axis_set);
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_refresh (int cycle, double time, double dt, int axis_set)
{
  // Update cycle,time,dt from Simulation

   // (should be updated already?)
  // -1 test due to p_refresh_face() calling p_refresh with axis_set != axis_all
  if (cycle != -1) {
    cycle_ = cycle;
    time_  = time;
    dt_    = dt;
  }

  refresh(axis_set);
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::refresh (int axis_set)
{

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

  Hierarchy *        hierarchy   = simulation->hierarchy();
  Boundary *         boundary    = simulation->boundary();
  const FieldDescr * field_descr = simulation->field_descr();

  double lower[3];
  hierarchy->lower(&lower[0],&lower[1],&lower[2]);
  double upper[3];
  hierarchy->upper(&upper[0],&upper[1],&upper[2]);

  bool ax = ((axis_set == axis_all) || (axis_set == axis_x)) && nx > 1;
  bool ay = ((axis_set == axis_all) || (axis_set == axis_y)) && ny > 1;
  bool az = ((axis_set == axis_all) || (axis_set == axis_z)) && nz > 1;

  if ( ax ) {
    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY
    boundary_face[axis_x][face_lower] = fabs(lower_[0]-lower[0]) < 1e-7;
    boundary_face[axis_x][face_upper] = fabs(upper_[0]-upper[0]) < 1e-7;
    if ( boundary_face[axis_x][face_lower] ) 
      boundary->enforce(field_descr,this,face_lower,axis_x);
    if ( boundary_face[axis_x][face_upper] ) 
      boundary->enforce(field_descr,this,face_upper,axis_x);
  }
  if ( ay ) {
    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY
    boundary_face[axis_y][face_lower] = fabs(lower_[1]-lower[1]) < 1e-7;
    boundary_face[axis_y][face_upper] = fabs(upper_[1]-upper[1]) < 1e-7;
    if ( boundary_face[axis_y][face_lower] ) 
      boundary->enforce(field_descr,this,face_lower,axis_y);
    if ( boundary_face[axis_y][face_upper] ) 
      boundary->enforce(field_descr,this,face_upper,axis_y);
  }
  if ( az ) {
    // COMPARISON INACCURATE FOR VERY SMALL BLOCKS NEAR BOUNDARY
    boundary_face[axis_z][face_lower] = fabs(lower_[2]-lower[2]) < 1e-7;
    boundary_face[axis_z][face_upper] = fabs(upper_[2]-upper[2]) < 1e-7;
    if ( boundary_face[axis_z][face_lower] ) 
      boundary->enforce(field_descr,this,face_lower,axis_z);
    if ( boundary_face[axis_z][face_upper] ) 
      boundary->enforce(field_descr,this,face_upper,axis_z);
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

  bool update_full = simulation->temp_update_full();

  if ( ax ) {
    // xp <<< xm
    if ( ! boundary_face[axis_x][face_lower] || periodic ) {
      field_face.load (field_descr, field_block(), axis_x, face_lower, 
		       update_full,update_full);
      block_array(ixm,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_x, face_upper, axis_set);

    }
    // xp >>> xm
    if ( ! boundary_face[axis_x][face_upper] || periodic ) {
      field_face.load (field_descr, field_block(), axis_x, face_upper, 
		       update_full,update_full);
      block_array(ixp,iy,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_x, face_lower, axis_set);
    }
  }
  if ( ay ) {
    // yp <<< ym
    if ( ! boundary_face[axis_y][face_lower] || periodic ) {
      field_face.load (field_descr, field_block(), axis_y, face_lower, 
		       update_full,update_full);
      block_array(ix,iym,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_y, face_upper, axis_set);
    }
    // yp >>> ym
    if ( ! boundary_face[axis_y][face_upper] || periodic ) {
      field_face.load (field_descr, field_block(), axis_y, face_upper, 
		       update_full,update_full);
      block_array(ix,iyp,iz).p_refresh_face 
	(field_face.size(), field_face.array(), axis_y, face_lower, axis_set);
    }
  }
  if ( az ) {
    // zp <<< zm
    if ( ! boundary_face[axis_z][face_lower] || periodic ) {
      field_face.load (field_descr, field_block(), axis_z, face_lower, 
		       update_full,update_full);
      block_array(ix,iy,izm).p_refresh_face 
	(field_face.size(), field_face.array(), axis_z, face_upper, axis_set);
    }
    // zp >>> zm
    if ( ! boundary_face[axis_z][face_upper] || periodic ) {
      field_face.load (field_descr, field_block(), axis_z, face_upper, 
		       update_full,update_full);
      block_array(ix,iy,izp).p_refresh_face 
	(field_face.size(), field_face.array(), axis_z, face_lower, axis_set);
    }
  }

  // NOTE: p_refresh_face() calls compute, but if no incoming faces
  // it will never get called.  So every block also calls
  // p_refresh_face() itself with a null array

  p_refresh_face (0,0,0,0, axis_set);

}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::p_refresh_face (int n, char * buffer,
			    int axis, int face, int axis_set)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  if ( n != 0) {

    // n == 0 is the "default" call from self to ensure p_refresh_face()
    // gets called at least once--required for subsequent compute() call

    FieldFace field_face(n, buffer);

    bool update_full = simulation->temp_update_full();

    field_face.store(simulation->field_descr(),
		     field_block(), 
		     axis_enum(axis), 
		     face_enum(face), 
		     update_full,
		     update_full);
  }

  //--------------------------------------------------
  // Count incoming faces
  //--------------------------------------------------

  Hierarchy * hierarchy = simulation->hierarchy();

  double lower[3];
  hierarchy->lower(&lower[0],&lower[1],&lower[2]);
  double upper[3];
  hierarchy->upper(&upper[0],&upper[1],&upper[2]);
  
  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  bool ax = (axis_set == axis_all || axis_set == axis_x) && nx > 1;
  bool ay = (axis_set == axis_all || axis_set == axis_y) && ny > 1;
  bool az = (axis_set == axis_all || axis_set == axis_z) && nz > 1;

  // Determine expected number of incoming Faces
  // (should be function, and possibly precomputed and stored since constant )

  int count = 1;  // self

  if (ax ) count +=2;
  if (ay ) count +=2;
  if (az ) count +=2;

  bool periodic = simulation->boundary()->is_periodic();

  if (! periodic && ax ) {
    if (fabs(lower_[0]-lower[0]) < 1e-7) --count;
    if (fabs(upper_[0]-upper[0]) < 1e-7) --count;
  }
  if (! periodic && ay ) {
    if (fabs(lower_[1]-lower[1]) < 1e-7) --count;
    if (fabs(upper_[1]-upper[1]) < 1e-7) --count;
  }
  if (! periodic && az ) {
    if (fabs(lower_[2]-lower[2]) < 1e-7) --count;
    if (fabs(upper_[2]-upper[2]) < 1e-7) --count;
  }

  //--------------------------------------------------
  // Compute
  //--------------------------------------------------

  if (axis_set == axis_all) {
    if (++count_refresh_face_ >= count) {
      count_refresh_face_ = 0;
#ifdef ORIGINAL_REFRESH  
      compute(axis_set);
#else
      prepare();
#endif
    }
  } else {
    switch (axis_set) {
    case axis_x:
      if (++count_refresh_face_x_ >= count) {
	count_refresh_face_x_ = 0;
	p_refresh (-1,0,0,axis_y);
      }
      break;
    case axis_y:
      if (++count_refresh_face_y_ >= count) {
	count_refresh_face_y_ = 0;
	p_refresh (-1,0,0,axis_z);
      }
      break;
    case axis_z:
      if (++count_refresh_face_z_ >= count) {
	count_refresh_face_z_ = 0;
#ifdef ORIGINAL_REFRESH  
	compute(axis_set); // axis_set ignored--used when compute() calls refresh()
#else
	prepare();
#endif
      }
      break;
    }
  }
}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

// SEE Simulation/simulation_charm_output.cpp for Block::p_write(int)

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::compute(int axis_set)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

#ifdef CONFIG_USE_PROJECTIONS
  double time_start = CmiWallTimer();
#endif

  // FieldDescr * field_descr = simulation->field_descr();
  // char buffer[10];
  // sprintf (buffer,"%03d-A",cycle_);
  // field_block()->print(field_descr,buffer,lower_,upper_);
  // field_block()->image(field_descr,"A",cycle_,
  // 			 thisIndex.x,thisIndex.y,thisIndex.z);

  for (size_t i = 0; i < simulation->num_method(); i++) {
    simulation->method(i) -> compute_block (this,time_,dt_);
  }

   // sprintf (buffer,"%03d-B",cycle_);
   // field_block()->print(field_descr,buffer,lower_,upper_);
   // field_block()->image(field_descr,"B",cycle_,
   // 		       thisIndex.x,thisIndex.y,thisIndex.z);

#ifdef CONFIG_USE_PROJECTIONS
  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  // Update Block cycle and time

  time_ += dt_;
  ++ cycle_ ;

  // prepare for next cycle: Timestep, Stopping, Monitor, Output

#ifdef ORIGINAL_REFRESH  
  prepare();
#else
  refresh(axis_set);
#endif

}
#endif /* CONFIG_USE_CHARM */

//======================================================================

  void Block::copy_(const Block & block) throw()
{

  // Create a copy of field_block_
  field_block_.resize(block.field_block_.size());
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (*(block.field_block_[i]));
  }
#ifdef CONFIG_USE_CHARM
  num_field_blocks_ = block.num_field_blocks_;
#endif
}

//----------------------------------------------------------------------
