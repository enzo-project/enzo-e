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
  count_coarsen_(0),
  count_adapt_(count_adapt),
  adapt_(adapt_unknown),
#ifdef CONFIG_USE_CHARM
  loop_refresh_(),
#endif
  face_level_()
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

  // Initialize neighbor

  face_level_.resize(27);
  for (int i=0; i<27; i++) face_level_[i] = -1;

  int na3[3];
  size_forest(&na3[0],&na3[1],&na3[2]);

  int v3[3];

  int ixp = (rank >= 1) ? 1 : 0;
  int iyp = (rank >= 2) ? 1 : 0;
  int izp = (rank >= 3) ? 1 : 0;
    
  int refresh_rank = this->simulation()->config()->field_refresh_rank;

  int icx=0,icy=0,icz=0;
  if (level_ > 0) {
    index_.child(level_,&icx,&icy,&icz);
  }

  bool periodic = this->simulation()->problem()->boundary()->is_periodic();
  for (int ix=-ixp; ix<=ixp; ix++) {
    bool isx = ((ix==0) || (ix==-1&&icx==1) || (ix==1&&icx==0));
    for (int iy=-iyp; iy<=iyp; iy++) {
      bool isy = ((iy==0) || (iy==-1&&icy==1) || (iy==1&&icy==0));
      for (int iz=-izp; iz<=izp; iz++) {
	bool isz = ((iz==0) || (iz==-1&&icz==1) || (iz==1&&icz==0));

	int face_rank = rank - (abs(ix)+abs(iy)+abs(iz));

	bool has_neighbor = ((refresh_rank <= face_rank && face_rank < rank) &&
			     ((level_ == 0) || (isx && isy && isz)));
	if (has_neighbor) {
	  Index index = index_.index_neighbor(ix,iy,iz,na3,periodic);
	  index.values(v3);
#ifdef CELLO_TRACE
	  index_.print("1 calling p_set_neighbor A");
	  index.print ("1 calling p_set_neighbor B");
#endif
	  int in3[3] = {ix,iy,iz};
	  TRACE("Calling p_set_neighbor");
	  p_set_neighbor (v3,in3);
	}
      }
    }
  }

  // 

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

//----------------------------------------------------------------------

void CommBlock::p_set_neighbor(const int v3[3], int in3[3])
{
  face_level_[IN3(in3)] = level_;
}

//----------------------------------------------------------------------

void CommBlock::p_delete_neighbor(const int v3[3],int in3[3])
{
  face_level_[IN3(in3)] = level_ - 1;
}

//----------------------------------------------------------------------

bool CommBlock::is_neighbor (const Index & index, int in3[3]) const
{ 
  return (face_level_[IN3(in3)] == level_);
}

//----------------------------------------------------------------------

void CommBlock::p_set_nibling(const int v3[3], int in3[3])
{
  face_level_[IN3(in3)] = level_ + 1;
}

//----------------------------------------------------------------------

void CommBlock::p_delete_nibling(const int v3[3], int in3[3])
{
  face_level_[IN3(in3)] = level_;
}

//----------------------------------------------------------------------

bool CommBlock::is_nibling (const Index & index, int in3[3]) const
{ 
  return (face_level_[IN3(in3)] == level_ + 1);
}

//======================================================================

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
