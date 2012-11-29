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
#ifdef CONFIG_USE_CHARM
 CProxy_Patch proxy_patch,
#endif
 int patch_id,
 int patch_rank,
 int num_field_blocks,
 CommBlock * comm_block
) throw ()
  :
#ifdef CONFIG_USE_CHARM
     count_refresh_face_(0),
     proxy_patch_(proxy_patch),
#endif
     num_field_blocks_(num_field_blocks),
     field_block_(),
     patch_id_(patch_id),
     patch_rank_(patch_rank),
     cycle_(0),
     time_(0),
     dt_(0),
     comm_block_(comm_block)
{ 
  DEBUG1("ID = %d",patch_id);
  DEBUG1("IP = %d",patch_rank);

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
#ifdef CONFIG_USE_CHARM
 CProxy_Patch proxy_patch,
#endif
 int patch_id,
 int patch_rank,
 int num_field_blocks,
 CommBlock * comm_block
 ) throw ()
  : count_refresh_face_(0),
    proxy_patch_(proxy_patch),
    num_field_blocks_(num_field_blocks),
    field_block_(),
    patch_id_(patch_id),
    patch_rank_(patch_rank),
    cycle_(0),
    time_(0),
    dt_(0),
    comm_block_(comm_block)
{ 

  DEBUG1("ID = %d",patch_id_);
  DEBUG1("IP = %d",patch_rank_);

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

//======================================================================
// CHARM FUNCTIONS
//======================================================================

#ifdef CONFIG_USE_CHARM
extern CProxy_SimulationCharm  proxy_simulation;
#endif /* CONFIG_USE_CHARM */

//======================================================================

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
  DEBUG("Block::prepare()");

  dt_block = timestep->evaluate(field_descr,this);
  DEBUG("Block::prepare()");

  // Reduce timestep to coincide with scheduled output if needed

  int index_output=0;
  while (Output * output = problem->output(index_output++)) {
    Schedule * schedule = output->schedule();
    dt_block = schedule->update_timestep(time_,dt_block);
  }

  DEBUG("Block::prepare()");

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

  CkCallback callback (CkIndex_Block::p_output(NULL), thisProxy);
  TRACE1("Calling contribute %d",2*sizeof(double));
  contribute( 2*sizeof(double), min_reduce, CkReduction::min_double, callback);

}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::refresh ()
{
  TRACE ("Block::refresh()");

  bool is_boundary[3][2];

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  Boundary * boundary = simulation->problem()->boundary();
  FieldDescr * field_descr = simulation->field_descr();
  
  bool periodic = boundary->is_periodic();

  CProxy_Block block_array = thisProxy;

  //--------------------------------------------------
  // Refresh
  //--------------------------------------------------

  int ix = thisIndex.x;
  int iy = thisIndex.y;
  int iz = thisIndex.z;

  int nbx = size_[0];
  int nby = size_[1];
  int nbz = size_[2];
  
  bool fx3[3],fy3[3],fz3[3];
  determine_boundary_
    (is_boundary,&fx3[0],&fx3[2],&fy3[0],&fy3[2],&fz3[0],&fz3[2]);
  fx3[1]=true;
  fy3[1]=true;
  fz3[1]=true;
  fx3[0] = fx3[0] && (periodic || ! is_boundary[axis_x][face_lower]);
  fx3[2] = fx3[2] && (periodic || ! is_boundary[axis_x][face_upper]);
  fy3[0] = fy3[0] && (periodic || ! is_boundary[axis_y][face_lower]);
  fy3[2] = fy3[2] && (periodic || ! is_boundary[axis_y][face_upper]);
  fz3[0] = fz3[0] && (periodic || ! is_boundary[axis_z][face_lower]);
  fz3[2] = fz3[2] && (periodic || ! is_boundary[axis_z][face_upper]);
  int ix3[3],iy3[3],iz3[3];
  ix3[0] = (ix - 1 + nbx) % nbx;
  iy3[0] = (iy - 1 + nby) % nby;
  iz3[0] = (iz - 1 + nbz) % nbz;
  ix3[1] = ix;
  iy3[1] = iy;
  iz3[1] = iz;
  ix3[2] = (ix + 1) % nbx;
  iy3[2] = (iy + 1) % nby;
  iz3[2] = (iz + 1) % nbz;
  
  // Refresh face ghost zones

  bool gx,gy,gz;
  gx = false;
  gy = false;
  gz = false;

  int fxl = 1;
  int fyl = (nby==1 && ! periodic) ? 0 : 1;
  int fzl = (nbz==1 && ! periodic) ? 0 : 1;

  for (int fx=-fxl; fx<=fxl; fx++) {
    for (int fy=-fyl; fy<=fyl; fy++) {
      for (int fz=-fzl; fz<=fzl; fz++) {
	int sum = abs(fx)+abs(fy)+abs(fz);
	if ((fx3[fx+1] && fy3[fy+1] && fz3[fz+1]) &&
	    ((sum==1 && field_descr->refresh_face(2)) ||
	     (sum==2 && field_descr->refresh_face(1)) ||
	     (sum==3 && field_descr->refresh_face(0)))) {

	  FieldFace field_face (field_block(),field_descr);

	  field_face.set_face(fx,fy,fz);
	  field_face.set_ghost(gx,gy,gz);
	  
	  DEBUG9("index %d %d %d  %d %d %d  %d %d %d",
		 index_[0],index_[1],index_[2],
		 ix3[fx+1],iy3[fy+1],iz3[fz+1],
		 fx,fy,fz);

	  int n; 
	  char * array;
	  field_face.load(&n, &array);

	  block_array(ix3[fx+1],iy3[fy+1],iz3[fz+1]).x_refresh (n, array, -fx,-fy,-fz);
	}
      }
    }
  }

  // NOTE: x_refresh() calls compute, but if no incoming faces
  // it will never get called.  So every block also calls
  // x_refresh() itself with a null array

  x_refresh (0,0,0, 0, 0);

}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::determine_boundary_
(
 bool is_boundary[3][2],
 bool * fxm,
 bool * fxp,
 bool * fym,
 bool * fyp,
 bool * fzm,
 bool * fzp
 )
{
  
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  Hierarchy * hierarchy   = simulation->hierarchy();

  double lower_h[3], upper_h[3];
  hierarchy->lower(&lower_h[0],&lower_h[1],&lower_h[2]);
  hierarchy->upper(&upper_h[0],&upper_h[1],&upper_h[2]);

  // return is_boundary[] array of faces on domain boundary

  is_on_boundary (lower_h,upper_h,is_boundary);

  int nx,ny,nz;
  field_block()->size (&nx,&ny,&nz);

  // Determine in which directions we need to communicate or update boundary

  if (fxm) *fxm = nx > 1;
  if (fxp) *fxp = nx > 1;
  if (fym) *fym = ny > 1;
  if (fyp) *fyp = ny > 1;
  if (fzm) *fzm = nz > 1;
  if (fzp) *fzp = nz > 1;
}

#endif

//----------------------------------------------------------------------


#ifdef CONFIG_USE_CHARM

void Block::update_boundary_
(
 bool is_boundary[3][2],
 bool fxm,
 bool fxp,
 bool fym,
 bool fyp,
 bool fzm,
 bool fzp
)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  Boundary * boundary = simulation->problem()->boundary();
  const FieldDescr * field_descr = simulation->field_descr();


  // Update boundaries

  if ( fxm && is_boundary[axis_x][face_lower] ) {
    boundary->enforce(field_descr,this,face_lower,axis_x);
  }
  if ( fxp && is_boundary[axis_x][face_upper] ) {
    boundary->enforce(field_descr,this,face_upper,axis_x);
  }
  if ( fym && is_boundary[axis_y][face_lower] ) {
    boundary->enforce(field_descr,this,face_lower,axis_y);
  }
  if ( fyp && is_boundary[axis_y][face_upper] ) {
    boundary->enforce(field_descr,this,face_upper,axis_y);
  }
  if ( fzm && is_boundary[axis_z][face_lower] ) {
    boundary->enforce(field_descr,this,face_lower,axis_z);
  }
  if ( fzp && is_boundary[axis_z][face_upper] ) {
    boundary->enforce(field_descr,this,face_upper,axis_z);
  }
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::compute()
{
  DEBUG ("Block::compute()");

  Simulation * simulation = proxy_simulation.ckLocalBranch();

 #ifdef CONFIG_USE_PROJECTIONS
   double time_start = CmiWallTimer();
 #endif

  FieldDescr * field_descr = simulation->field_descr();

  int index_method = 0;
  while (Method * method = simulation->problem()->method(index_method++)) {
    method -> compute_block (field_descr,this);
  }

 #ifdef CONFIG_USE_PROJECTIONS
   traceUserBracketEvent(10,time_start, CmiWallTimer());
 #endif

  // Update Block cycle and time to Simulation time and cycle

  set_cycle (cycle_ + 1);

  set_time  (time_  + dt_);
  
  // prepare for next cycle: Timestep, Stopping, Monitor, Output

  DEBUG ("Block::compute() calling refresh()");
  refresh();

}
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Block::initial(Simulation * simulation)
{
  TRACE("Block::initial()");

  FieldDescr * field_descr = simulation->field_descr();

  // Initialize the block

  allocate(field_descr);

  // Set the Block cycle and time to match Simulation's

  TRACE("CommBlock::p_initial Setting time");
  set_cycle(simulation->cycle());
  set_time (simulation->time());
  set_dt   (simulation->dt());

  // Perform any additional initialization for derived class 

  initialize ();

  // Apply the initial conditions 

  Initial * initial = simulation->problem()->initial();

  initial->enforce_block(this,field_descr, simulation->hierarchy());

}
#endif

//----------------------------------------------------------------------

void Block::output(Simulation * simulation)
{
  DEBUG("Block::output()");

  double * min_reduce = (double * )msg->getData();

  double dt_patch   = min_reduce[0];
  bool   stop_patch = min_reduce[1] == 1.0 ? true : false;

  delete msg;

  set_dt   (dt_patch);

  // WARNING: assumes one patch

  simulation->update_state(cycle_,time_,dt_patch,stop_patch);
 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // ??? HOW IS cycle_ and time_ update on all processors ensured before index() calls
  // Simulation::p_output()?  Want last block?
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  

  // "root" block calls Simulation::p_output()
  //  if (index() == 0) {
  //    proxy_simulation.p_output();
  //  }

  // Wait for all blocks to check in before calling Simulation::p_output()
  // for next output

}

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

void Block::allocate (const FieldDescr * field_descr) throw()
{ 
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i]->allocate_array(field_descr,true);
    //    field_block_[i]->allocate_array(field_descr);
    //    field_block_[i]->allocate_ghosts(field_descr);
  }
}

//----------------------------------------------------------------------

