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
  face_level_(),
  face_level_new_(),
  child_face_level_(),
  child_face_level_new_(),
  count_coarsen_(0),
  level_count_(index.level()),
  adapt_step_(num_adapt_steps),
  adapt_(adapt_unknown),
  next_phase_(phase_stopping),
  coarsened_(false),
  delete_(false),
  is_leaf_(true)
{

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

  block_ = new Block  (nx, ny, nz, num_field_blocks,
		       xm, xp, ym, yp, zm, zp);
  block_->allocate(field_descr);

  child_block_ = NULL;

  // Update state

  set_state (cycle,time,dt,stop_);

  // Perform any additional initialization for derived class 

  initialize ();

  for (int i=0; i<3; i++) {
    count_sync_[i] = 0;
    max_sync_[i] = 0;
  }

  // Initialize neighbor face levels

  const int rank = simulation()->dimension();
  if (num_face_level == 0) {

    face_level_.resize(27);
    child_face_level_.resize(NC(rank)*27);

    for (int i=0; i<27; i++) face_level_[i] = 0;
    for (int i=0; i<NC(rank)*27; i++) child_face_level_[i] = 0;

  } else {

    face_level_.resize(num_face_level);
    child_face_level_.resize(NC(rank)*num_face_level);

    for (int i=0; i<num_face_level; i++) face_level_[i] = face_level[i];

  }

  //  face_level_new_.resize(face_level_.size());
  face_level_new_ = face_level_;
  //  child_face_level_new_.resize(child_face_level_.size());
  child_face_level_new_ = child_face_level_;

  initialize_child_face_levels_();

  debug_faces_("CommBlock");

  const int level = this->level();

  int na3[3];
  size_forest(&na3[0],&na3[1],&na3[2]);

  int icx=0,icy=0,icz=0;
  if (level > 0) index_.child(level,&icx,&icy,&icz);

  if (narray != 0) {
    
    // copy any input data
    FieldDescr * field_descr = simulation()->field_descr();
    FieldFace field_face (block()->field_block(),field_descr);

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

    // --------------------------------------------------
    // ENTRY: #1 CommBlock::CommBlock() -> CommBlock::q_adapt_exit()
    // ENTRY: quiescence if level > 0
    // --------------------------------------------------
    CkStartQD (CkCallback(CkIndex_CommBlock::q_adapt_exit(), 
			  thisProxy[thisIndex]));
    // --------------------------------------------------
  }
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
#ifdef TEMP_NEW_REFINE
  p | level_desired_;
#endif
  p | cycle_;
  p | time_;
  p | dt_;
  p | stop_;
  //  p | neighbor_index_;
  p | index_initial_;
  p | children_;
  p | loop_refresh_;
  p | sync_coarsen_;
  PUParray(p,count_sync_, PHASE_SYNC_SIZE);
  PUParray(p,max_sync_, PHASE_SYNC_SIZE);
  p | face_level_;
  p | face_level_new_;
  p | child_face_level_;
  p | child_face_level_new_;
  p | count_coarsen_;
  p | level_count_;
  p | adapt_step_;
  p | adapt_;
  p | next_phase_;
  p | coarsened_;
  p | delete_;

}

//----------------------------------------------------------------------

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

}

//----------------------------------------------------------------------

Simulation * CommBlock::simulation() const
{ return proxy_simulation.ckLocalBranch(); }

//----------------------------------------------------------------------

std::string CommBlock::name() const throw()
{
  int dim = simulation()->dimension();
  return std::string("Block-") + index_.bit_string(level(),dim);
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

//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

// void CommBlock::delete_child(Index index)
// {
//   for (size_t i=0; i<children_.size(); i++) {
//     // erase by replacing occurences with self
//     if (children_[i] == index) children_[i] = index_;
//   }
// }

//======================================================================

// int CommBlock::count_neighbors() const
// {
//   if (! is_leaf()) return 0;
//   int level = this->level();
//   const int rank         = simulation()->dimension();
//   const int rank_refresh = simulation()->config()->field_refresh_rank;
//   const bool periodic    = simulation()->problem()->boundary()->is_periodic();
//   ItFace it_face (rank,rank_refresh);
//   int num_neighbors = 0;
//   int ic3[3] = {0,0,0};
//   int of3[3];
//   while (it_face.next(of3)) {
//     Index index_neighbor = neighbor_(of3);
//     const int level_face = face_level (of3);
//     if (level_face == level) { // SAME
//       ++num_neighbors;
//     } else if (level_face == level - 1) { // COARSE
//       index_.child (level,&ic3[0],&ic3[1],&ic3[2]);
//       int op3[3];
//       parent_face_(op3,of3,ic3);
//       if (op3[0]==of3[0] && op3[1]==of3[1] && op3[2]==of3[2]) 
// 	++num_neighbors;
//     } else if (level_face == level + 1) { // FINE
//       const int if3[3] = {-of3[0],-of3[1],-of3[2]};
//       ItChild it_child(rank,if3);
//       while (it_child.next(ic3)) 
// 	++num_neighbors;
//     } else {
//       ERROR2 ("CommBlock::count_neighbors()",
// 	      "level_face %d level %d",
// 	      level_face,level);
//     }
//   }
// #ifdef CELLO_DEBUG
//   char buffer[255];
//   sprintf (buffer,"count_neighbors = %d",num_neighbors);
//   index_.print(buffer,-1,2,false,simulation());
// #endif
//   return num_neighbors;

// }
//----------------------------------------------------------------------

void CommBlock::loop_limits_refresh_(int ifacemin[3], int ifacemax[3])
  const throw()
{

  Boundary * boundary = simulation()->problem()->boundary();

  // which faces need to be refreshed?
  bool on_boundary[3][2];
  is_on_boundary (on_boundary);

  bool periodic = boundary->is_periodic();
  if (periodic) {
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
	on_boundary[axis][face] = false;
      }
    }
  }

  // set face loop limits accordingly
  ifacemin[0] = on_boundary[0][0] ? 0 : -1;
  ifacemin[1] = on_boundary[1][0] ? 0 : -1;
  ifacemin[2] = on_boundary[2][0] ? 0 : -1;
  ifacemax[0] = on_boundary[0][1] ? 0 : 1;
  ifacemax[1] = on_boundary[1][1] ? 0 : 1;
  ifacemax[2] = on_boundary[2][1] ? 0 : 1;

  int rank = simulation()->dimension();
  if (rank < 2) ifacemin[1] = ifacemax[1] = 0;
  if (rank < 3) ifacemin[2] = ifacemax[2] = 0;

}

//----------------------------------------------------------------------

void CommBlock::loop_limits_nibling_ 
( int ic3m[3],int ic3p[3], const int if3[3]) const throw()
{
  int rank = simulation()->dimension();

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

Index CommBlock::neighbor_ 
(
 const int of3[3],
 Index *   ind
 ) const
{
  Index index = (ind != 0) ? (*ind) : index_;

  int na3[3];
  size_forest (&na3[0],&na3[1],&na3[2]);
  const bool periodic  = simulation()->problem()->boundary()->is_periodic();
  Index in = index.index_neighbor (of3[0],of3[1],of3[2],na3,periodic);
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

      index_.print(mesg,-1,2,false,simulation());

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
	fprintf (fp_debug,(ic3[1]==1) ? "%d " : "  ",face_level_new(if3));
#endif
	PARALLEL_PRINTF ((ic3[1]==1) ? "%d " : "  ",face_level_new(if3));
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
	    fprintf (fp_debug,"%d ",child_face_level_new(ic3,if3));
#endif
	    PARALLEL_PRINTF ("%d ",child_face_level_new(ic3,if3));
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

