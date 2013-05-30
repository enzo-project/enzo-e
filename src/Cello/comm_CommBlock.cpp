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
 int num_field_blocks,
 int count_adapt,
 bool initial,
 int cycle, double time, double dt,
 int narray, char * array, int op_array,
 int num_face_level, int * face_level,
 bool testing
 ) throw ()
  : 
#ifndef CONFIG_USE_CHARM
  simulation_(simulation),
#endif
  index_(index),
  index_initial_(0),
  level_(index_.level()),
  children_(),
#ifdef CONFIG_USE_CHARM
  loop_refresh_(),
#endif
  face_level_(),
  count_coarsen_(0),
  count_adapt_(count_adapt),
  adapt_(adapt_unknown)
{ 
#ifdef CELLO_TRACE
  index.print ("CommBlock::CommBlock");
  printf("CommBlock::CommBlock(n(%d %d %d)  num_field_blocks %d  count_adapt %d  initial %d)\n",
	 nx,ny,nz,num_field_blocks,count_adapt,initial);

  printf("CommBlock::CommBlock  n (%d %d %d)\n",nx,ny,nz);
  printf("CommBlock::CommBlock  l %d\n",level_);
#endif


  int ibx,iby,ibz;
  index.array(&ibx,&iby,&ibz);

  double xm,ym,zm;
  lower(&xm,&ym,&zm);
  double xp,yp,zp;
  upper(&xp,&yp,&zp);

  FieldDescr * field_descr = this->simulation()->field_descr();

  block_ = new Block  (nx, ny, nz, num_field_blocks,
		       xm, xp, ym, yp, zm, zp);
  // Allocate block data
  block_->allocate(field_descr);
  child_block_ = NULL;

#ifdef CONFIG_USE_CHARM

  // Call virtual functions to update state

  set_cycle(cycle);
  set_time (time);
  set_dt   (dt);

#endif

  // Perform any additional initialization for derived class 

  int rank = this->simulation()->dimension();

  initialize ();

  // Initialize neighbor face levels

  if (num_face_level == 0) {

    face_level_.resize(27);
    for (int i=0; i<27; i++) face_level_[i] = 0;
  } else {
    face_level_.resize(num_face_level);
    for (int i=0; i<num_face_level; i++) face_level_[i] = face_level[i];
  }

  int na3[3];
  size_forest(&na3[0],&na3[1],&na3[2]);

  int icx=0,icy=0,icz=0;
  if (level_ > 0) {
    index_.child(level_,&icx,&icy,&icz);
  }

  if (narray != 0) {
    // copy any input data
    FieldDescr * field_descr = this->simulation()->field_descr();
    FieldFace field_face (block()->field_block(),field_descr);

    //    set "face" to full FieldBlock
    field_face.set_face(0,0,0);

    //    set array operation if any
    switch (op_array) {
    case op_array_restrict:
      {
	Restrict * restrict = this->simulation()->problem()->restrict();
	field_face.set_restrict(restrict,icx,icy,icz);
      }
      break;
    case op_array_prolong:
      {
	Prolong * prolong = this->simulation()->problem()->prolong();
	field_face.set_prolong(prolong,icx,icy,icz);
      }
      break;
    default:
      break;
    }

    field_face.store(narray,array);
  }

#ifdef CONFIG_USE_CHARM

  if (! testing) {

    // Count CommBlocks on each processor
   
    ((SimulationCharm *)simulation())->insert_block();

  }

  if (initial) apply_initial_();

#endif /* CONFIG_USE_CHARM */

}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void CommBlock::pup(PUP::er &p)
{
  TRACEPUP;

  CBase_CommBlock::pup(p);

  p | *block_;
  p | *child_block_;
  p | index_;
  p | cycle_;
  p | time_;
  p | dt_;
  p | index_initial_;
  p | level_;
  p | children_;
  p | count_coarsen_;
  p | count_adapt_;
  p | adapt_;
  p | loop_refresh_;
  p | face_level_;

}

#endif /* CONFIG_USE_CHARM */

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

  TRACE("CommBlock::apply_initial_()");
  simulation()->performance()->start_region(perf_initial);
  FieldDescr * field_descr = simulation()->field_descr();

  // Apply initial conditions

  index_initial_ = 0;
  Problem * problem = simulation()->problem();
  while (Initial * initial = problem->initial(index_initial_++)) {
    initial->enforce_block(this,field_descr, simulation()->hierarchy());
  }
  simulation()->performance()->stop_region(perf_initial);
}
#endif

//----------------------------------------------------------------------

CommBlock::~CommBlock() throw ()
{ 
#ifdef CONFIG_USE_CHARM

  if (level_ > 0) {

    // Send restricted data to parent 

    // <duplicated code: refactor me!>
    FieldBlock * field_block = block()->field_block();
    FieldDescr * field_descr = simulation()->field_descr();
    FieldFace field_face (field_block,field_descr);

    //    set "face" to full FieldBlock
    field_face.set_face(0,0,0);

    //    set restriction
    Restrict * restrict = simulation()->problem()->restrict();
    int icx,icy,icz;
    index_.child(level_,&icx,&icy,&icz);
    field_face.set_restrict(restrict,icx,icy,icz);

    //    load array with data
    int n; 
    char * array;
    field_face.load(&n,&array);
    // </duplicated code>

    //    send data to parent

    thisProxy[index_.index_parent()].x_refresh_child(n,array,icx,icy,icz);
    
  }
  if (block_) delete block_;
  block_ = 0;
  if (child_block_) delete child_block_;
  child_block_ = 0;

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

//----------------------------------------------------------------------

bool CommBlock::is_child (const Index & index) const
{ 
  for (size_t i=0; i<children_.size(); i++) {
    if (children_[i] == index) return true;
  }
  return false;
}

//----------------------------------------------------------------------

void CommBlock::delete_child(Index index)
{
  for (size_t i=0; i<children_.size(); i++) {
    // erase by replacing occurences with self
    if (children_[i] == index) children_[i] = index_;
  }
}

//======================================================================

void CommBlock::loop_limits_refresh_
(int * ifxm, int * ifym, int * ifzm, int * ifxp, int * ifyp, int * ifzp)
  const throw()
{

  Boundary * boundary = simulation()->problem()->boundary();

  // which faces need to be refreshed?
  bool on_boundary[3][2];
  is_on_boundary (on_boundary);
  TRACE6("REFRESH is_on_boundary %d %d %d  %d %d %d",
	 on_boundary[0][0],on_boundary[1][0],on_boundary[2][0],
	 on_boundary[0][1],on_boundary[1][1],on_boundary[2][1]);

  bool periodic = boundary->is_periodic();
  if (periodic) {
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
	on_boundary[axis][face] = false;
      }
    }
  }

  // set face loop limits accordingly
  (*ifxm) = on_boundary[0][0] ? 0 : -1;
  (*ifym) = on_boundary[1][0] ? 0 : -1;
  (*ifzm) = on_boundary[2][0] ? 0 : -1;
  (*ifxp) = on_boundary[0][1] ? 0 : 1;
  (*ifyp) = on_boundary[1][1] ? 0 : 1;
  (*ifzp) = on_boundary[2][1] ? 0 : 1;

  int rank = simulation()->dimension();
  if (rank < 2) (*ifym) = (*ifyp) = 0;
  if (rank < 3) (*ifzm) = (*ifzp) = 0;

}

//----------------------------------------------------------------------

void CommBlock::loop_limits_nibling_ 
(int *icxm, int *icym, int *iczm,
 int *icxp, int *icyp, int *iczp,
 int ifx, int ify, int ifz) const throw()
{
  int rank = simulation()->dimension();

  (*icxm) = (ifx == 0) ? 0 : (ifx+1)/2;
  (*icxp) = (ifx == 0) ? 1 : (ifx+1)/2;
  (*icym) = (ify == 0) ? 0 : (ify+1)/2;
  (*icyp) = (ify == 0) ? 1 : (ify+1)/2;
  (*iczm) = (ifz == 0) ? 0 : (ifz+1)/2;
  (*iczp) = (ifz == 0) ? 1 : (ifz+1)/2;
  if (rank < 2) (*icym) = (*icyp) = 0;
  if (rank < 3) (*iczm) = (*iczp) = 0;
}

//----------------------------------------------------------------------

void CommBlock::copy_(const CommBlock & comm_block) throw()
{
  block_->copy_(*comm_block.block());
  if (child_block_) child_block_->copy_(*comm_block.child_block());

  cycle_ = comm_block.cycle_;
  time_  = comm_block.time_;
  dt_    = comm_block.dt_;
  level_ = comm_block.level_;
  count_adapt_ = comm_block.count_adapt_;
  adapt_ = comm_block.adapt_;
}

//----------------------------------------------------------------------

void CommBlock::is_on_boundary (bool is_boundary[3][2]) const throw()
{

  Boundary * boundary = simulation()->problem()->boundary();
  bool periodic = boundary->is_periodic();

  int n3[3];
  size_forest (&n3[0],&n3[1],&n3[2]);
  
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      is_boundary[axis][face] = 
	index_.is_on_boundary(axis,face,n3[axis],periodic);
    }
  }
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


//----------------------------------------------------------------------

void CommBlock::debug_faces_(const char * buffer, int * faces)
{
#ifdef CELLO_TRACE
  if (!faces) return;
  int imp[3] = {-1,1,0};
  int i0p[3] = { 0,1,0};
  int ipp[3] = { 1,1,0};
  int im0[3] = {-1,0,0};
  int i00[3] = { 0,0,0};
  int ip0[3] = { 1,0,0};
  int imm[3] = {-1,-1,0};
  int i0m[3] = { 0,-1,0};
  int ipm[3] = { 1,-1,0};
  index_.print(buffer);
  printf ("%s %2d %2d %2d\n",
   	  buffer, faces[IF3(imp)], faces[IF3(i0p)], faces[IF3(ipp)]);
  printf ("%s %2d %2d %2d\n",
   	  buffer, faces[IF3(im0)], faces[IF3(i00)], faces[IF3(ip0)]);
  printf ("%s %2d %2d %2d\n",
   	  buffer, faces[IF3(imm)], faces[IF3(i0m)], faces[IF3(ipm)]);
#endif
}
