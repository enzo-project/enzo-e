// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Feb 28 13:22:26 PST 2011
/// @brief    Implementation of the CommBlock object
/// @todo     Remove hierarchy dependency via Initial--only need domain bounds

#include "cello.hpp"

#include "mesh.hpp"
#include "main.hpp"
#include "charm_simulation.hpp"

#include "enzo.hpp" /* temp */

//----------------------------------------------------------------------

CommBlock::CommBlock
(
#ifndef CONFIG_USE_CHARM
 Simulation * simulation,
#endif
 Index index,
 int nx, int ny, int nz,             // Block cells
 int level,                          // Block level
 int num_field_blocks,
 int count_adapt,
 bool testing
 ) throw ()
  : 
#ifndef CONFIG_USE_CHARM
  simulation_(simulation),
#endif
  index_(index),
  cycle_(0),
  time_(0),
  dt_(0),
  index_initial_(0),
  level_(level),
  children_(),
  neighbors_(),
  niblings_(),
  count_adapt_(count_adapt)
{ 

  TRACE("ENTER CommBlock::CommBlock()");

#ifdef CONFIG_USE_CHARM
  index.print("create");
#endif

  TRACE3("CommBlock::CommBlock  n (%d %d %d)",nx,ny,nz);
  TRACE1("CommBlock::CommBlock  l %d",level);

  int ibx,iby,ibz;
  index.array(&ibx,&iby,&ibz);

  double xm,ym,zm;
  lower(&xm,&ym,&zm);
  double xp,yp,zp;
  upper(&xp,&yp,&zp);

  block_ = new Block  (nx, ny, nz, num_field_blocks,
		       xm, xp, ym, yp, zm, zp);

#ifdef CONFIG_USE_CHARM
  FieldDescr * field_descr = simulation()->field_descr();

  // Set the CommBlock cycle and time to match Simulation's

  TRACE("CommBlock::p_initial Setting time");
  set_cycle(simulation()->cycle());
  set_time (simulation()->time());
  set_dt   (simulation()->dt());
  // Allocate block data

  block()->allocate(field_descr);

#endif
  // Perform any additional initialization for derived class 

  int rank = this->simulation()->dimension();

  initialize ();

  if (level == 0) {
    WARNING("CommBlock::CommBlock",
	    "Ignoring periodic versus non-periodic b.c.");
    int na3[3];
    size_forest(&na3[0],&na3[1],&na3[2]);
    for (int axis = 0; axis < rank; axis++) {
      for (int face = 0; face < 2; face ++) {
	p_set_neighbor (index_.index_neighbor(axis,face,na3[axis]));
	
      }
    }
  }
#ifdef CONFIG_USE_CHARM

   if (! testing) {
     // Count CommBlocks on each processor
   
     TRACE1 ("simulation = %p",simulation());
     TRACE1 ("proxy_simulation = %p",&proxy_simulation);
     ((SimulationCharm *)simulation())->insert_block();
   }

  sync_refresh_.stop() = count_refresh_();

  apply_initial_();

#endif /* CONFIG_USE_CHARM */

  //

}

//----------------------------------------------------------------------

void CommBlock::initialize_
(
 int nx, int ny, int nz,
 bool testing
 )
 {
 }

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
void CommBlock::apply_initial_() throw ()
{

  FieldDescr * field_descr = simulation()->field_descr();

  // Apply initial conditions

  index_initial_ = 0;
  Problem * problem = simulation()->problem();
  while (Initial * initial = problem->initial(index_initial_++)) {
    initial->enforce_block(this,field_descr, simulation()->hierarchy());
  }
}
#endif

//----------------------------------------------------------------------

CommBlock::~CommBlock() throw ()
{ 
#ifdef CONFIG_USE_CHARM

  if (block_) delete block_;
  block_ = 0;

  ((SimulationCharm *)simulation())->delete_block();
#endif

}

//----------------------------------------------------------------------

CommBlock::CommBlock(const CommBlock & block) throw ()
/// @param     block  Object being copied
{
  copy_(block);
}

//----------------------------------------------------------------------

CommBlock & CommBlock::operator = (const CommBlock & block) throw ()
/// @param     block  Source object of the assignment
/// @return    The target assigned object
{
  copy_(block);
  return *this;
}

//----------------------------------------------------------------------

void CommBlock::index_forest (int * ix, int * iy, int * iz) const throw ()
{
  index_.array(ix,iy,iz);
}

//----------------------------------------------------------------------

void CommBlock::size_forest (int * nx=0, int * ny=0, int * nz=0) const throw ()
{
  simulation()->hierarchy()->num_blocks(nx,ny,nz);

  TRACE3 ("size_forest = %d %d %d",*nx,*ny,*nz);
}

//----------------------------------------------------------------------

void CommBlock::lower
(double * xm, double * ym, double * zm) const throw ()
{
  int  ix, iy, iz;
  int  nx, ny, nz;

  index_global (&ix,&iy,&iz,&nx,&ny,&nz);

  double xdm, ydm, zdm;
  double xdp, ydp, zdp;
  simulation()->hierarchy()->lower(&xdm,&ydm,&zdm);
  simulation()->hierarchy()->upper(&xdp,&ydp,&zdp);

  double ax = 1.0*ix/nx;
  double ay = 1.0*iy/ny;
  double az = 1.0*iz/nz;

  double xbm = (1.0-ax)*xdm + ax*xdp;
  double ybm = (1.0-ay)*ydm + ay*ydp;
  double zbm = (1.0-az)*zdm + az*zdp;

  if (xm) (*xm) = xbm;
  if (ym) (*ym) = ybm;
  if (zm) (*zm) = zbm;
  TRACE6 ("DEBUG LOWER %d %d %d  %f %f %f",ix,iy,iz,*xm,*ym,*zm);
}

//----------------------------------------------------------------------

void CommBlock::upper
(double * xp, double * yp, double * zp) const throw ()
{
  int  ix, iy, iz;
  int  nx, ny, nz;

  index_global (&ix,&iy,&iz,&nx,&ny,&nz);

  double xdm, ydm, zdm;
  double xdp, ydp, zdp;
  simulation()->hierarchy()->lower(&xdm,&ydm,&zdm);
  simulation()->hierarchy()->upper(&xdp,&ydp,&zdp);

  double ax = 1.0*(ix+1)/nx;
  double ay = 1.0*(iy+1)/ny;
  double az = 1.0*(iz+1)/nz;

  double xbp = (1.0-ax)*xdm + ax*xdp;
  double ybp = (1.0-ay)*ydm + ay*ydp;
  double zbp = (1.0-az)*zdm + az*zdp;

  if (xp) (*xp) = xbp;
  if (yp) (*yp) = ybp;
  if (zp) (*zp) = zbp;
}

//----------------------------------------------------------------------

void CommBlock::index_global
( int *ix, int *iy, int *iz,
  int *nx, int *ny, int *nz ) const
{
  index_forest(ix,iy,iz);
  size_forest (nx,ny,nz);
  TRACE6("DEBUG INDEX A %d %d %d  %d %d %d",
	 *ix,*iy,*iz,*nx,*ny,*nz);

  Index index = this->index();

  int level = index.level();
  for (int i=0; i<level; i++) {
    int bx,by,bz;
    index.child(i+1,&bx,&by,&bz);
    if (ix) (*ix) = ((*ix) << 1) | bx;
    if (iy) (*iy) = ((*iy) << 1) | by;
    if (iz) (*iz) = ((*iz) << 1) | bz;
    if (nx) (*nx) <<= 1;
    if (ny) (*ny) <<= 1;
    if (nz) (*nz) <<= 1;
  }
  TRACE6("DEBUG INDEX B %d %d %d  %d %d %d",
	 *ix,*iy,*iz,*nx,*ny,*nz);
}

//======================================================================
// MPI FUNCTIONS
//======================================================================

#ifdef CONFIG_USE_CHARM
#else /* CONFIG_USE_CHARM */

void CommBlock::refresh_ghosts(const FieldDescr * field_descr,
			       const Hierarchy * hierarchy,
			       int fx, int fy, int fz,
			       int index_field_set) throw()
{
  int ibx,iby,ibz;

  index_forest(&ibx,&iby,&ibz);

  block_->field_block(index_field_set)
    -> refresh_ghosts (field_descr,
		       hierarchy->group_process(),
		       hierarchy->layout(),
		       ibx,iby,ibz, fx,fy,fz);
}

#endif /* CONFIG_USE_CHARM */

#ifdef CONFIG_USE_CHARM

void CommBlock::determine_boundary_
(
 bool is_boundary[3][2],
 bool * fxm, bool * fxp,
 bool * fym, bool * fyp,
 bool * fzm, bool * fzp
 )
{
  
  // return is_boundary[] array of faces on domain boundary

  is_on_boundary (is_boundary);

  int nx,ny,nz;
  block_->field_block()->size (&nx,&ny,&nz);

  // Determine in which directions we need to communicate or update boundary

  if (fxm) *fxm = (nx > 1);
  if (fxp) *fxp = (nx > 1);
  if (fym) *fym = (ny > 1);
  if (fyp) *fyp = (ny > 1);
  if (fzm) *fzm = (nz > 1);
  if (fzp) *fzp = (nz > 1);
}

#endif

//----------------------------------------------------------------------


#ifdef CONFIG_USE_CHARM

void CommBlock::update_boundary_ ()
{

  bool is_boundary[3][2];
  bool fxm,fxp,fym,fyp,fzm,fzp;

  determine_boundary_(is_boundary,&fxm,&fxp,&fym,&fyp,&fzm,&fzp);

  Boundary * boundary = simulation()->problem()->boundary();
  const FieldDescr * field_descr = simulation()->field_descr();


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

//======================================================================

void CommBlock::copy_(const CommBlock & comm_block) throw()
{
  block_->copy_(*comm_block.block());

  cycle_ = comm_block.cycle_;
  time_  = comm_block.time_;
  dt_    = comm_block.dt_;
  level_ = comm_block.level_;
  count_adapt_ = comm_block.count_adapt_;

#ifdef CONFIG_USE_CHARM
  sync_refresh_ = comm_block.sync_refresh_;
#endif
  
}

//----------------------------------------------------------------------

void CommBlock::is_on_boundary (bool is_boundary[3][2]) throw()
{
  int ix,iy,iz;
  index_forest(&ix,&iy,&iz);
  int nx,ny,nz;
  size_forest(&nx,&ny,&nz);

  is_boundary[0][0] = (ix == 0);
  is_boundary[1][0] = (iy == 0);
  is_boundary[2][0] = (iz == 0);
  is_boundary[0][1] = (ix == (nx - 1));
  is_boundary[1][1] = (iy == (ny - 1));
  is_boundary[2][1] = (iz == (nz - 1));
}

//----------------------------------------------------------------------

Simulation * CommBlock::simulation() const
{
#ifdef CONFIG_USE_CHARM
  return proxy_simulation.ckLocalBranch();
#else /* CONFIG_USE_CHARM */
  return simulation_;
#endif /* CONFIG_USE_CHARM */
}

