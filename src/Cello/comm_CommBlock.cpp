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

// #define DEBUG_ADAPT

//----------------------------------------------------------------------

CommBlock::CommBlock
(
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int num_adapt_steps,
 int cycle, double time, double dt,
 int narray, char * array, int op_array,
 int num_face_level, int * face_level,
 bool testing
 ) throw ()
  :
  index_(index),
  stop_(0),
  index_initial_(0),
  children_(),
  loop_refresh_(),
  sync_coarsen_(),
  count_sync_(),
  max_sync_(),
  face_level_curr_(),
  face_level_next_(),
  child_face_level_curr_(),
  child_face_level_next_(),
  count_coarsen_(0),
  adapt_step_(num_adapt_steps),
  adapt_(adapt_unknown),
  next_phase_(phase_stopping),
  coarsened_(false),
  delete_(false),
  is_leaf_(true),
  age_(0),
  face_level_last_(),
  name_(name()),
  method_(0)
{

  // Enable Charm++ AtSync() dynamic load balancing
  usesAtSync = CmiTrue;

#ifdef CELLO_DEBUG
  index_.print("CommBlock()",-1,2,false,simulation());
#endif

  int ibx,iby,ibz;
  index.array(&ibx,&iby,&ibz);

  double xm,ym,zm;
  lower(&xm,&ym,&zm);
  double xp,yp,zp;
  upper(&xp,&yp,&zp);

  FieldDescr * field_descr = simulation()->field_descr();

  // Allocate block data

  block_ = new Block  (field_descr,
		       nx, ny, nz, 
		       num_field_blocks,
		       xm,xp, ym,yp, zm,zp);


  block_->allocate(field_descr);

  child_block_ = NULL;

  // Update state

  set_state (cycle,time,dt,stop_);

  // Perform any additional initialization for derived class 

  initialize ();

  const int rank = this->rank();
  
  sync_coarsen_.set_stop(NC(rank));
  sync_coarsen_.clear();

  for (int i=0; i<3; i++) {
    count_sync_[i] = 0;
    max_sync_[i] = 0;
  }

  // Initialize neighbor face levels

  face_level_last_.resize(27*8);

  if (num_face_level == 0) {

    face_level_curr_.resize(27);
    child_face_level_curr_.resize(NC(rank)*27);

    for (int i=0; i<27; i++) face_level_curr_[i] = 0;

  } else {

    face_level_curr_.resize(num_face_level);
    child_face_level_curr_.resize(NC(rank)*num_face_level);

    for (int i=0; i<num_face_level; i++) face_level_curr_[i] = face_level[i];

  }

  for (int i=0; i<face_level_last_.size(); i++) face_level_last_[i] = 0;
  for (int i=0; i<child_face_level_curr_.size(); i++) child_face_level_curr_[i] = 0;

  initialize_child_face_levels_();

  face_level_next_ = face_level_curr_;
  child_face_level_next_ = child_face_level_curr_;

  const int level = this->level();

  int na3[3];
  size_forest(&na3[0],&na3[1],&na3[2]);

  int icx=0,icy=0,icz=0;
  if (level > 0) index_.child(level,&icx,&icy,&icz);

  if (narray != 0) {
    
    // copy any input data
    FieldFace field_face (block()->field_block());

    //    set "face" to full FieldBlock
    field_face.set_face(0,0,0);
    field_face.set_ghost(true,true,true);

    //    set array operation if any

    Problem * problem = simulation()->problem();

    switch (op_array) {

    case op_array_restrict:
      field_face.set_restrict(problem->restrict(),icx,icy,icz);  
      break;

    case op_array_prolong:
      field_face.set_prolong(problem->prolong(),icx,icy,icz);
      break;

    default:
      break;

    }

    field_face.store(narray,array);

  }

  if (! testing) ((SimulationCharm *)simulation())->insert_block();

  int initial_cycle = simulation()->config()->initial_cycle;
  bool is_first_cycle = (initial_cycle == cycle);

  if (is_first_cycle) {
    apply_initial_();
  } else if (level > 0) {

    thisProxy.doneInserting();
    control_sync (sync_adapt_end,"quiescence");

  }

  debug_faces_("CommBlock()");

}

//----------------------------------------------------------------------

void CommBlock::pup(PUP::er &p)
{
  TRACEPUP;

  CBase_CommBlock::pup(p);

  bool up = p.isUnpacking();

  if (up) block_ = new Block;
  p | *block_;

  // child_block_ may be NULL
  bool allocated=(child_block_ != NULL);
  p|allocated;
  if (allocated) {
    if (up) child_block_=new Block;
    p|*child_block_;
  } else {
    child_block_ = NULL;
  }

  p | index_;
  p | level_next_;
  p | cycle_;
  p | time_;
  p | dt_;
  p | stop_;
  p | index_initial_;
  p | children_;
  p | loop_refresh_;
  p | sync_coarsen_;
  PUParray(p,count_sync_, SYNC_SIZE);
  PUParray(p,max_sync_,   SYNC_SIZE);
  p | face_level_curr_;
  p | face_level_next_;
  p | child_face_level_curr_;
  p | child_face_level_next_;
  p | count_coarsen_;
  p | adapt_step_;
  p | adapt_;
  p | next_phase_;
  p | coarsened_;
  p | delete_;
  p | is_leaf_;
  p | age_;
  p | face_level_last_;
  p | name_;

  // SKIP method_: initialized when needed

  if (up) debug_faces_("PUP");
}

//----------------------------------------------------------------------

const FieldDescr * CommBlock::field_descr() throw()
{ return block_->field_block()->field_descr(); }

//----------------------------------------------------------------------

ItFace CommBlock::it_face(const int * ic3,
			  const int * if3) throw()
{
  const Config * config = simulation()->config();
  int rank = this->rank();
  int refresh_rank = config->field_refresh_rank;
  int n3[3];
  size_forest(&n3[0],&n3[1],&n3[2]);
  bool periodic[3][2];
  periodicity(periodic);
  return ItFace (rank,refresh_rank,periodic,n3,index_,ic3,if3);
}

//======================================================================

void CommBlock::apply_initial_() throw ()
{

  TRACE("CommBlock::apply_initial_()");

  performance_switch_(perf_initial,__FILE__,__LINE__);

  FieldDescr * field_descr = simulation()->field_descr();

  // Apply initial conditions

  index_initial_ = 0;
  Problem * problem = simulation()->problem();
  while (Initial * initial = problem->initial(index_initial_++)) {
    initial->enforce_block(this,field_descr, simulation()->hierarchy());
  }
  //  performance_stop_(perf_initial);
}

//----------------------------------------------------------------------

CommBlock::~CommBlock() throw ()
{ 
#ifdef CELLO_DEBUG
  index_.print("~CommBlock()",-1,2,false,simulation());
#endif

  const int level = this->level();

  if (level > 0) {

    // Send restricted data to parent 

    int ichild[3];
    index_.child(level,ichild,ichild+1,ichild+2);

    int n; 
    char * array;
    int iface[3]={0,0,0};
    bool lghost[3]={true,true,true};

    FieldFace * field_face = 
      load_face_(&n,&array,iface,ichild,lghost,op_array_restrict);

    const Index index_parent = index_.index_parent();

    // --------------------------------------------------
    // ENTRY: #2 CommBlock::~CommBlock()-> CommBlock::x_refresh_child()
    // ENTRY: parent if level > 0
    // --------------------------------------------------
    thisProxy[index_parent].x_refresh_child(n,array,ichild);
    // --------------------------------------------------

    delete field_face;
  }
  if (block_) delete block_;
  block_ = 0;
  if (child_block_) delete child_block_;
  child_block_ = 0;

  ((SimulationCharm *)simulation())->delete_block();

  thisProxy.doneInserting();

}

//----------------------------------------------------------------------

CommBlock::CommBlock (CkMigrateMessage *m) : CBase_CommBlock(m)
{ 
  ((SimulationCharm *)simulation())->insert_block();
};

//----------------------------------------------------------------------

Simulation * CommBlock::simulation() const
{ return proxy_simulation.ckLocalBranch(); }

//----------------------------------------------------------------------

int CommBlock::rank() const
{ return simulation()->rank(); }

//----------------------------------------------------------------------

std::string CommBlock::name() const throw()
{
  const int rank = this->rank();
  int nb3[3] = {1,1,1};
  simulation()->hierarchy()->blocking(nb3,nb3+1,nb3+2);
  int nb  = std::max( std::max (nb3[0],nb3[1]),nb3[2]);
  int bits = 0;
  while (nb/=2) ++bits;
  return "B" + index_.bit_string(level(),rank,bits);
}

//----------------------------------------------------------------------

void CommBlock::size_forest (int * nx, int * ny, int * nz) const throw ()
{  simulation()->hierarchy()->num_blocks(nx,ny,nz); }

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

void CommBlock::cell_width 
(double * dx, double * dy, double * dz) const throw()
{ 
  double xm,ym,zm;
  lower(&xm,&ym,&zm);
  double xp,yp,zp;
  upper(&xp,&yp,&zp);
  block()->field_block()->cell_width(xm,xp,dx, ym,yp,dy, zm,zp,dz);
}

//----------------------------------------------------------------------

void CommBlock::index_global
( int *ix, int *iy, int *iz,
  int *nx, int *ny, int *nz ) const
{
  
  index_forest(ix,iy,iz);
  size_forest (nx,ny,nz);

  Index index = this->index();

  const int level = this->level();

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

//----------------------------------------------------------------------

void CommBlock::determine_boundary_
(
 bool bndry[3][2],
 bool * fxm, bool * fxp,
 bool * fym, bool * fyp,
 bool * fzm, bool * fzp
 )
{
  // return bndry[] array of faces on domain boundary

  is_on_boundary (bndry);

  int nx,ny,nz;
  block()->field_block()->size (&nx,&ny,&nz);

  // Determine in which directions we need to communicate or update boundary

  if (fxm && bndry[0][0]) *fxm = (nx > 1);
  if (fxp && bndry[0][1]) *fxp = (nx > 1);
  if (fym && bndry[1][0]) *fym = (ny > 1);
  if (fyp && bndry[1][1]) *fyp = (ny > 1);
  if (fzm && bndry[2][0]) *fzm = (nz > 1);
  if (fzp && bndry[2][1]) *fzp = (nz > 1);
}

//----------------------------------------------------------------------

void CommBlock::update_boundary_ ()
{
  bool is_boundary[3][2];
  bool fxm=0,fxp=0,fym=0,fyp=0,fzm=0,fzp=0;

  determine_boundary_(is_boundary,&fxm,&fxp,&fym,&fyp,&fzm,&fzp);

  const FieldDescr * field_descr = simulation()->field_descr();

  int index = 0;
  Problem * problem = simulation()->problem();
  Boundary * boundary;

  while ((boundary = problem->boundary(index++))) {
    // Update boundaries
    if ( fxm ) boundary->enforce(field_descr,this,face_lower,axis_x);
    if ( fxp ) boundary->enforce(field_descr,this,face_upper,axis_x);
    if ( fym ) boundary->enforce(field_descr,this,face_lower,axis_y);
    if ( fyp ) boundary->enforce(field_descr,this,face_upper,axis_y);
    if ( fzm ) boundary->enforce(field_descr,this,face_lower,axis_z);
    if ( fzp ) boundary->enforce(field_descr,this,face_upper,axis_z);
  }
}

//----------------------------------------------------------------------

void CommBlock::periodicity (bool p32[3][2]) const
{
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      p32[axis][face] = false;
    }
  }
  int index_boundary = 0;
  Problem * problem = simulation()->problem();
  Boundary * boundary;
  while ( (boundary = problem->boundary(index_boundary++)) ) {
    boundary->periodicity(p32);
  }
}

//----------------------------------------------------------------------

void CommBlock::loop_limits_nibling_ 
( int ic3m[3],int ic3p[3], const int if3[3]) const throw()
{
  const int rank = this->rank();

  ic3m[0] = (if3[0] == 0) ? 0 : (if3[0]+1)/2;
  ic3p[0] = (if3[0] == 0) ? 1 : (if3[0]+1)/2;
  ic3m[1] = (if3[1] == 0) ? 0 : (if3[1]+1)/2;
  ic3p[1] = (if3[1] == 0) ? 1 : (if3[1]+1)/2;
  ic3m[2] = (if3[2] == 0) ? 0 : (if3[2]+1)/2;
  ic3p[2] = (if3[2] == 0) ? 1 : (if3[2]+1)/2;
  if (rank < 2) ic3m[1] = ic3p[1] = 0;
  if (rank < 3) ic3m[2] = ic3p[2] = 0;
}

//----------------------------------------------------------------------

void CommBlock::facing_child_(int jc3[3], const int ic3[3], const int if3[3]) const
{
  jc3[0] = if3[0] ? 1 - ic3[0] : ic3[0];
  jc3[1] = if3[1] ? 1 - ic3[1] : ic3[1];
  jc3[2] = if3[2] ? 1 - ic3[2] : ic3[2];
  // if (if3[0]) jc3[0] = 1 - jc3[0];
  // if (if3[0]) jc3[0] = 1 - jc3[0];
  // if (if3[0]) jc3[0] = 1 - jc3[0];
  // if (if3[1]==-1) jc3[1] = 1;
  // if (if3[2]==-1) jc3[2] = 1;
  // if (if3[0]==+1) jc3[0] = 0;
  // if (if3[1]==+1) jc3[1] = 0;
  // if (if3[2]==+1) jc3[2] = 0;
  TRACE9("facing_child %d %d %d  child %d %d %d  face %d %d %d",
	 jc3[0],jc3[1],jc3[2],ic3[0],ic3[1],ic3[2],if3[0],if3[1],if3[2]);
}

//----------------------------------------------------------------------

void CommBlock::copy_(const CommBlock & comm_block) throw()
{
  block_->copy_(*comm_block.block());
  if (child_block_) child_block_->copy_(*comm_block.child_block());

  cycle_      = comm_block.cycle_;
  time_       = comm_block.time_;
  dt_         = comm_block.dt_;
  stop_       = comm_block.stop_;
  adapt_step_ = comm_block.adapt_step_;
  adapt_      = comm_block.adapt_;
  next_phase_ = comm_block.next_phase_;
  coarsened_  = comm_block.coarsened_;
  delete_     = comm_block.delete_;
}

//----------------------------------------------------------------------

void CommBlock::is_on_boundary (bool is_boundary[3][2]) const throw()
{

  // Boundary * boundary = simulation()->problem()->boundary();
  // bool periodic = boundary->is_periodic();

  int n3[3];
  size_forest (&n3[0],&n3[1],&n3[2]);
  
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      is_boundary[axis][face] = 
	index_.is_on_boundary(axis,face,n3[axis]);
    }
  }
}

//----------------------------------------------------------------------

Index CommBlock::neighbor_ 
(
 const int of3[3],
 Index *   ind
 ) const
{
  Index index = (ind != 0) ? (*ind) : index_;

  int na3[3];
  size_forest (&na3[0],&na3[1],&na3[2]);
  // const bool periodic  = simulation()->problem()->boundary()->is_periodic();
  Index in = index.index_neighbor (of3,na3);
  return in;
}

//----------------------------------------------------------------------

void CommBlock::performance_start_
(int index_region, std::string file, int line)
{
  simulation()->performance()->start_region(index_region,file,line);
}

//----------------------------------------------------------------------

void CommBlock::performance_stop_
(int index_region, std::string file, int line)
{
  simulation()->performance()->stop_region(index_region,file,line);
}

//----------------------------------------------------------------------

void CommBlock::performance_switch_
(int index_region, std::string file, int line)
{
  simulation()->performance()->switch_region(index_region,file,line);
}

//----------------------------------------------------------------------

void CommBlock::debug_faces_(const char * mesg)
{
#ifndef DEBUG_ADAPT
  return;
#endif

#ifdef CELLO_DEBUG
  FILE * fp_debug = simulation()->fp_debug();
#endif
  TRACE_ADAPT(mesg);
  int if3[3] = {0};
  int ic3[3] = {0};

  for (ic3[1]=1; ic3[1]>=0; ic3[1]--) {
    for (if3[1]=1; if3[1]>=-1; if3[1]--) {

      index_.print(mesg,-1,2,true,simulation());

      for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	fprintf (fp_debug,(ic3[1]==1) ? "%d " : "  ",face_level(if3));
#endif
	PARALLEL_PRINTF ((ic3[1]==1) ? "%d " : "  ",face_level(if3));
      }
#ifdef CELLO_DEBUG
      fprintf (fp_debug,"| ");
#endif
      PARALLEL_PRINTF ("| ") ;
      for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	fprintf (fp_debug,(ic3[1]==1) ? "%d " : "  ",face_level_next(if3));
#endif
	PARALLEL_PRINTF ((ic3[1]==1) ? "%d " : "  ",face_level_next(if3));
      }
#ifdef CELLO_DEBUG
      fprintf (fp_debug,"| ");
#endif
      PARALLEL_PRINTF ("| ");
      for (ic3[0]=0; ic3[0]<2; ic3[0]++) {
	for (if3[0]=-1; if3[0]<=1; if3[0]++) {
	  for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	    fprintf (fp_debug,"%d ",child_face_level(ic3,if3));
#endif
	    PARALLEL_PRINTF ("%d ",child_face_level(ic3,if3));
	  }
	}
      }
#ifdef CELLO_DEBUG
      fprintf (fp_debug,"| ");
#endif
      PARALLEL_PRINTF ("| ");
      for (ic3[0]=0; ic3[0]<2; ic3[0]++) {
	for (if3[0]=-1; if3[0]<=1; if3[0]++) {
	  for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	    fprintf (fp_debug,"%d ",child_face_level_next(ic3,if3));
#endif
	    PARALLEL_PRINTF ("%d ",child_face_level_next(ic3,if3));
	  }
	}
      }
#ifdef CELLO_DEBUG
      fprintf (fp_debug,"\n");
      fflush(fp_debug);
#endif
      PARALLEL_PRINTF ("\n");
      fflush(stdout);

    }
  }
}

