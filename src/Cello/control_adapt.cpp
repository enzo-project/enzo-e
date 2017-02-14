// See LICENSE_CELLO file for license and copyright information

/// @file     control_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions.  
/// @ingroup  Control
///
/// This file controls adaptive mesh refinement on a distributed
/// forest of octrees.

//--------------------------------------------------
// #define DEBUG_ADAPT
// #define DEBUG_NEW_REFRESH
//--------------------------------------------------

#ifdef DEBUG_ADAPT

#define CELLO_TRACE

#   define PUT_LEVEL(INDEX_SEND,INDEX_RECV,IC3,IF3,LEVEL_NOW,LEVEL_NEW,MSG) \
  {									\
    char buffer [256];							\
    int nb3[3] = {2,2,2};						\
    CkPrintf ("%s %s -> B%s"						\
	     " [%d => %d] if3 %2d %2d %2d  ic3 %d %d %d\n",		\
	      name().c_str(),MSG,INDEX_RECV.bit_string(INDEX_RECV.level(),rank(),nb3).c_str(),LEVEL_NOW,LEVEL_NEW, \
	     IF3[0],IF3[1],IF3[2],IC3[0],IC3[1],IC3[2]);		\
    check_child_(IC3,"PUT_LEVEL",__FILE__,__LINE__);			\
    check_face_(IF3,"PUT_LEVEL",__FILE__,__LINE__);			\
    thisProxy[INDEX_RECV].p_adapt_recv_level				\
      (INDEX_SEND,IC3,IF3,LEVEL_NOW,LEVEL_NEW);				\
  }
#else /* DEBUG_ADAPT */
#   define PUT_LEVEL(INDEX_SEND,INDEX_RECV,IC3,IF3,LEVEL_NOW,LEVEL_NEW,MSG) \
  {									\
    thisProxy[INDEX_RECV].p_adapt_recv_level				\
      (INDEX_SEND,IC3,IF3,LEVEL_NOW,LEVEL_NEW);				\
  }
#endif /* DEBUG_ADAPT */

#ifdef DEBUG_ADAPT

#   ifdef CELLO_TRACE
#      define trace(A) CkPrintf ("%s:%d %s TRACE %s\n",			\
				 __FILE__,__LINE__,name_.c_str(),A); fflush(stdout)
#   else
#      define trace(A) /*  NULL */
#   endif

#else

#   define trace(A) /*  NULL */

#endif

//--------------------------------------------------

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//======================================================================

/// @brief First function in the adapt phase: apply local refinement criteria.
/// 
/// adapt_begin_() computes the desired refinement level for the block
/// using adapt_compute_desired_level_(), after which it calls
/// adapt_called_() with nearest-neighbor synchronization.

void Block::adapt_begin_()

{
  trace("adapt_begin 1");

  simulation()->set_phase(phase_adapt);

  int level_maximum = simulation()->config()->mesh_max_level;

  level_next_ = adapt_compute_desired_level_(level_maximum);

  control_sync (CkIndex_Block::p_adapt_called(),sync_neighbor,0);

}

//----------------------------------------------------------------------

/// @brief Second step of the adapt phase: tell neighbors desired level.
///
/// Call adapt_send_level() to send neighbors desired
/// levels, after which adapt_next_() is called with quiescence
/// detection.
void Block::adapt_called_()
{
  trace("adapt_called 2");

  adapt_send_level();

  control_sync (CkIndex_Main::p_adapt_next(),sync_quiescence);
}

//----------------------------------------------------------------------

/// @brief Third step of the adapt phase: coarsen or refine according
/// to desired level.
///
/// Call update_levels_() to finalize face and child face levels,
/// then, if a leaf, refine or coarsen according to desired level
/// determined in adapt_called_().  Afterward, all Blocks call
/// adapt_end_().
void Block::adapt_next_()
{
  debug_faces_("adapt_next");
  trace("adapt_next 3");

#ifdef DEBUG_ADAPT
  {
    char buffer[255];
    if (level() != level_next_) {
      CkPrintf ("%s is leaf %d level %d -> %d\n",
		name().c_str(),is_leaf(),level(),level_next_);
    }
  }
#endif

  update_levels_();

  if (is_leaf()) {
    if (level() < level_next_) adapt_refine_();
    if (level() > level_next_) adapt_coarsen_();
  }

  control_sync (CkIndex_Main::p_adapt_end(),sync_quiescence);
}

//----------------------------------------------------------------------

/// @brief Fourth step of the adapt phase: delete self if Block has
/// been coarsened
///
/// This step deletes itself if it has been coarsened in this adapt
/// phase, then exits the adapt phase by directly calling adapt_exit_().
/// This is a separate phase since the quiescence call of this function
/// from the previous adapt_next_() step includes Block's that have
/// been deleted.  
void Block::adapt_end_()
{
  trace("adapt_end 4");

  if (index_.is_root()) {
    thisProxy.doneInserting();
  }

  if (delete_) {
#ifdef DEBUG_ADAPT
  CkPrintf ("%s DESTROY\n",name().c_str());
  fflush(stdout);
#endif
    ckDestroy();
    return;
  }

  for (size_t i=0; i<face_level_last_.size(); i++)
    face_level_last_[i] = -1;

  const int rank = this->rank();
  sync_coarsen_.set_stop(NUM_CHILDREN(rank));
  sync_coarsen_.reset();

  const int initial_cycle = simulation()->config()->initial_cycle;
  const bool is_first_cycle = (initial_cycle == cycle());
  const int level_maximum = simulation()->config()->mesh_max_level;

  bool adapt_again = (is_first_cycle && (adapt_step_++ < level_maximum));

  if (adapt_again) {
    control_sync (CkIndex_Main::p_adapt_enter(),sync_quiescence);
  } else {
    control_sync (CkIndex_Main::p_adapt_exit(),sync_quiescence);
  }

}

//----------------------------------------------------------------------

/// @brief Return whether the adapt phase should be called this cycle.
bool Block::do_adapt_()
{
  int adapt_interval = simulation()->config()->adapt_interval;

  return ((adapt_interval && ((cycle_ % adapt_interval) == 0)));
}

//----------------------------------------------------------------------

/// @brief Determine whether this Block should refine, coarsen, or
/// stay the same.
///
/// Return if not a leaf; otherwise, apply all Refine refinement
/// criteria to the Block, and set level_desired accordingly:
/// level+1 if it needs to refine, level - 1 if it can coarsen,
/// or level.
///
/// @param[in]  level_maximum   Maximum level to refine
///
/// @return The desired level based on local refinement criteria.
int Block::adapt_compute_desired_level_(int level_maximum)
{
  if (! is_leaf()) return adapt_same;

  adapt_ = adapt_unknown;

  int level = this->level();
  int level_desired = level;

  Problem * problem = simulation()->problem();
  Refine * refine;

  int index_refine = 0;
  while ((refine = problem->refine(index_refine++))) {

    adapt_ = std::max(adapt_,refine->apply(this));

  }

  const int initial_cycle = simulation()->config()->initial_cycle;
  const bool is_first_cycle = (initial_cycle == cycle());

  if (adapt_ == adapt_coarsen && level > 0 && ! is_first_cycle) 
    level_desired = level - 1;
  else if (adapt_ == adapt_refine  && level < level_maximum) 
    level_desired = level + 1;
  else {
    adapt_ = adapt_same;
    level_desired = level;
  }

  return level_desired;
}

//----------------------------------------------------------------------

void Block::adapt_refine_()
{
  Monitor * monitor = simulation()->monitor();
  if (monitor->is_verbose()) {
    char buffer [80];
    int v3[3];
    index().values(v3);
    sprintf (buffer,"Block %s (%x %x %x) is refining",name().c_str(),
	     v3[0],v3[1],v3[2]);
    monitor->print("Adapt",buffer);
  }

#ifdef DEBUG_ADAPT
  CkPrintf ("%s REFINE\n",name().c_str());
  fflush(stdout);
#endif

  adapt_ = adapt_unknown;

  const int rank = this->rank();
  
  int nx,ny,nz;
  data()->field_data()->size(&nx,&ny,&nz);

  std::vector<int> field_list;

  // First scatter particles to children first to avoid multiple passes
  
  ParticleData * particle_list[8] = {0};

  ParticleDescr * p_descr = simulation()->particle_descr();

  const int nc = NUM_CHILDREN(rank);

  for (int i=0; i<nc; i++) {
    particle_list[i] = new ParticleData;
    particle_list[i]->allocate(p_descr);
  }

  Particle particle = data()->particle();

  // partition particles into 2x2x2 array to send to children

  particle_scatter_children_ (particle_list,particle);
  
  // For each new child

  ItChild it_child (rank);
  int ic3[3];
  while (it_child.next(ic3)) {

    Index index_child = index_.index_child(ic3);

    // If child doesn't exist yet                                               

    if ( ! is_child_(index_child) ) {

      // Create FieldFace for interpolating field data to child ic3[]

      int if3[3] = {0,0,0};
      bool lg3[3] = {true,true,true};
      
      FieldFace * field_face = create_face 
	(if3,ic3,lg3, refresh_fine,field_list);

      DataMsg * data_msg = NULL;

      int narray = 0;
      char * array = 0;
      int num_field_data = 1;

      // Create data message object to send

      data_msg = new DataMsg;

      data_msg -> set_field_face (field_face,false);
      data_msg -> set_field_data (data()->field_data(),false);
      ParticleData * p_data = new ParticleData(*particle_list[IC3(ic3)]);
      data_msg -> set_particle_data (p_data,true);

      const Factory * factory = simulation()->factory();

      // Create the child object with interpolated data

      factory->create_block 
	(
	 data_msg,
	 &thisProxy, index_child,
	 nx,ny,nz,
	 num_field_data,
	 adapt_step_,
	 cycle_,time_,dt_,
	 narray, array, refresh_fine,
	 27,&child_face_level_curr_[27*IC3(ic3)],
	 simulation());

      delete [] array;
      array = 0;

      children_.push_back(index_child);

    }
  }

  for (int i=0; i<nc; i++) {
    delete particle_list[i];
  }

  // Delete particles since relocated to children

  int count = 0;
  int nt = particle.num_types();
  for (int it=0; it<nt; it++) {
    int nb = particle.num_batches(it);
    for (int ib=0; ib<nb; ib++) {
      count += particle.delete_particles (it, ib);
    }
  }
  simulation()->monitor_delete_particles(count);
  
  is_leaf_ = false;
#ifdef DEBUG_ADAPT
  CkPrintf ("%s adapt_refine is_leaf <- 0\n",name().c_str());
  fflush(stdout);
#endif
}

//----------------------------------------------------------------------

void Block::particle_scatter_children_ (ParticleData * particle_list[],
					Particle particle)
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("DEBUG_NEW_REFRESH particle_scatter_children\n");
#endif

  const int rank = this->rank();
  const int npa = NUM_CHILDREN(rank);

  // get Block bounds 
  double xm,ym,zm;
  double xp,yp,zp;
  lower(&xm,&ym,&zm);
  upper(&xp,&yp,&zp);

  // find block center (x0,y0,z0)
  const double x0 = 0.5*(xm+xp);
  const double y0 = 0.5*(ym+yp);
  const double z0 = 0.5*(zm+zp);

  // for each particle type to be moved

  const int nt = particle.num_types();

  int count = 0;
  for (int it=0; it<nt; it++) {

#ifdef DEBUG_NEW_REFRESH
    CkPrintf ("DEBUG_NEW_REFRESH scatter type %d\n",it);
#endif
    const int ia_x  = particle.attribute_position(it,0);

    // (...positions may use absolute coordinates (float) or
    // block-local coordinates (int))
    const bool is_float = 
      (cello::type_is_float(particle.attribute_type(it,ia_x)));

    // (...stride may be != 1 if particle attributes are interleaved)
    const int d  = particle.stride(it,ia_x);

    // ...for each batch of particles

    const int nb = particle.num_batches(it);

    for (int ib=0; ib<nb; ib++) {

      const int np = particle.num_particles(it,ib);

      // ...all particles will be moved
      bool mask[np];
      for (int ip=0; ip<np; ip++) {
	mask[ip] = true;
      }

      // ...extract particle position arrays

      double xa[np],ya[np],za[np];
      particle.position(it,ib,xa,ya,za);

      // ...initialize mask used for scatter and delete
      // ...and corresponding particle indices

      int index[np];

      if (is_float) {

	// absolute coordinates

	for (int ip=0; ip<np; ip++) {

#ifdef DEBUG_NEW_REFRESH
    CkPrintf ("DEBUG_NEW_REFRESH scatter particle %d\n",ip);
#endif

	  double x = xa[ip*d];
	  double y = ya[ip*d];
	  double z = za[ip*d];

	  int ix = (rank >= 1) ? ( (x < x0) ? 0 : 1) : 0;
	  int iy = (rank >= 2) ? ( (y < y0) ? 0 : 1) : 0;
	  int iz = (rank >= 3) ? ( (z < z0) ? 0 : 1) : 0;

	  // save index of ip'th particle's destination Particle object 
	  index[ip] = ix + 2*(iy + 2*iz);
	}

      } else {
	ERROR("Block::particle_scatter_children_",
	      "Relative (integer) positions not supported");
      }

      // ...scatter particles to particle array
      particle.scatter (it,ib, np, mask, index, npa, particle_list);
      // ... delete scattered particles
      count += particle.delete_particles (it,ib,mask);
    }
  }
  simulation()->monitor_delete_particles(count);
}

//----------------------------------------------------------------------

void Block::adapt_delete_child_(Index index_child)
{
#ifdef DEBUG_ADAPT
  int nb3[3] = {2,2,2};
  CkPrintf ("%s deleting child %s\n",
	    name().c_str(), index_child.bit_string(index_child.level(),rank(),nb3).c_str());
  fflush(stdout);
#endif
  thisProxy[index_child].p_adapt_delete();

  if (sync_coarsen_.next()) {
    children_.clear();
  }
}

//----------------------------------------------------------------------

void Block::adapt_send_level()
{
  if (!is_leaf()) return;

  const int level = this->level();
  const int min_face_rank = 
    simulation()->config()->adapt_min_face_rank;
  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);
  int of3[3];

  while (it_neighbor.next()) {
    Index index_neighbor = it_neighbor.index();
    int ic3[3];
    it_neighbor.child(ic3);
    it_neighbor.face(of3);
    PUT_LEVEL (index_,index_neighbor,ic3,of3,level,level_next_,"send");
  }
}

//----------------------------------------------------------------------

/// @brief Entry function for receiving desired level of a neighbor
///
/// @param index_send      mesh index of the calling neighbor
/// @param ic3             child indices of neighbor if it's in a finer level
/// @param if3             face (inward) shared with neighbor
/// @param level_face_curr neighbor's current level
/// @param level_face_new  neighbor's desired level
///
/// level_face_curr
/// level_face_new
/// level_next
/// level_next_
/// level
void Block::p_adapt_recv_level
(
 Index index_send,
 int ic3[3],
 int if3[3], 
 int level_face_curr,
 int level_face_new
 )
{
  performance_start_(perf_adapt_update);
  
  if (index_send.level() != level_face_curr) {
    PARALLEL_PRINTF 
      ("%d level mismatch between index_send (%d) and level_face_curr (%d)",
       __LINE__,index_send.level(), level_face_curr);
    index_.print("index_",-1,2,false,simulation());
    index_send.print("index_",-1,2,false,simulation());
  }

  // note face_level_last_ initialized as -1, in which case won't skip
  const bool skip_face_update = 
    (level_face_new <= face_level_last_[ICF3(ic3,if3)]);

#ifdef DEBUG_ADAPT
  {
    char buffer [255];
    int nb3[3] = {2,2,2};				
    CkPrintf ("%s %s <- B%s"
	      " [%d => %d] if3 %2d %2d %2d  ic3 %d %d %d [%d] %s\n",
	      name().c_str(),"recv",
	      index_send.bit_string(index_send.level(),rank(),nb3).c_str(),
	      level_face_curr,level_face_new,
	      if3[0],if3[1],if3[2],				
	      ic3[0],ic3[1],ic3[2], face_level_last_[ICF3(ic3,if3)],
	      skip_face_update ? "SKIP" : "");
    fflush(stdout);
  }
#endif

  if (skip_face_update) {
    performance_stop_(perf_adapt_update);
    performance_start_(perf_adapt_update_sync);
    return;
  }

  face_level_last_[ICF3(ic3,if3)] = level_face_new;

  int level_next = level_next_;

  const int level = this->level();

  const int of3[3] = {-if3[0],-if3[1],-if3[2] };

  if ( ! is_leaf()) {

    ERROR1 ("Block::p_adapt_recv_level()",
	    "Block %s is not a leaf",
	    name().c_str());
  }


  if (level_face_curr == level) {

    adapt_recv(of3,ic3,level_face_new,0);

  } else if ( level_face_curr == level + 1 ) {

    adapt_recv(of3,ic3,level_face_new,+1);

  } else if ( level_face_curr == level - 1 ) {

    adapt_recv(of3,ic3,level_face_new,-1);

  } else  {

    WARNING2 ("Block::notify_neighbor()",
	      "level %d and face level %d differ by more than 1",
	      level,level_face_curr);
  }

  // If this block wants to coarsen, then 
  //
  //    1. all siblings must be able to coarsen as well, and
  //
  //    2. all non-sibling (nephew) must not be going to a level
  //       finer than the current level (otherwise a level jump will
  //       occur).
  //
  // If either of these cases is true, then change the desired level
  // to the current level (neither coarsen nor refine) and re-send
  // desired level to neighbors

  const bool is_coarsening = (level_next < level);

  // The calling block is a sibling if it has the same parent

  bool is_sibling = false;

  if (level > 0 && index_send.level() > 0) {
    Index parent      = index_.    index_parent();
    Index send_parent = index_send.index_parent();
    is_sibling = (parent == send_parent);
  }

  // The calling block is a nephew if it is a child of a sibling

  bool is_nephew = false;

  if (level > 0 && index_send.level() > 1) {
    Index parent             = index_    .index_parent(); 
    Index send_parent_parent = index_send.index_parent().index_parent();
    is_nephew = (parent == send_parent_parent);
  }

  const bool is_finer_neighbor = (level_face_new > level_next);

  // If want to coarsen, then ensure that all siblings want to coarsen too

  if (is_coarsening) {

    // assume we can coarsen
    bool can_coarsen = true;

    // cannot if neighbor is sibling and not coarsening
    if (is_sibling && is_finer_neighbor) {
      can_coarsen = false;
    }

    // cannot if sibling has children
    if (is_nephew) {
      can_coarsen = false;
    }

    if (! can_coarsen) {
      level_next = level;
    }
  }

  // restrict new level to within 1 of neighbor
  level_next = std::max(level_next,level_face_new - 1);
	  
  // notify neighbors if level_next has changed

  if (level_next != level_next_) {
    ASSERT2 ("Block::p_adapt_recv_level()",
	     "level_next %d level_next_ %d\n", level_next,level_next_,
	     level_next > level_next_);
    level_next_ = level_next;
    adapt_send_level();
  }
  performance_stop_(perf_adapt_update);
  performance_start_(perf_adapt_update_sync);
}

//----------------------------------------------------------------------

void Block::adapt_recv 
( const int of3[3], const int ic3[3], int level_face_new, int level_relative )
{

  const int rank = this->rank();
  const int min_face_rank = 
    simulation()->config()->adapt_min_face_rank;

  if (level_relative == 0 || level_relative == +1) {
    // RECV-SAME: Face and level are received from unique
    // neighbor.  Unique face level is updated, and levels on
    // possibly multiple faces of multiple children are updated.

    set_face_level_next (of3, level_face_new);

    ItChild it_child (rank,of3);
    int jc3[3];

    while (it_child.next(jc3)) {

      int kf3[3];

      Index index_child = index_.index_child(jc3);
      ItFace it_face = this->it_face(min_face_rank,index_child,jc3,of3);

      if (level_relative == 0) {
	while (it_face.next()) {
	  it_face.face(kf3);
	  set_child_face_level_next(jc3,kf3,level_face_new);
	} 
      } else if (level_relative == +1) {
	while (it_face.next()) {

	  it_face.face(kf3);

	  Index in = neighbor_(kf3,&index_child);

	  Index index_neighbor = neighbor_(of3).index_child(ic3);

	  if (in == index_neighbor) {
	    set_child_face_level_next(jc3,kf3,level_face_new);
	  }
	}
      }
    }

  } else if (level_relative == -1) {

    // RECV-COARSE: Face and level are received from unique
    // neighbor.  Possibly multiple faces of block are updated
    // corresponding to the coarse neighbor's face.  Levels of
    // possibly multiple faces of possibly multiple child faces are
    // updated.
    ItFace it_face = this->it_face(min_face_rank,index_,ic3,of3);

    int jf3[3];
    while (it_face.next()) {
      it_face.face(jf3);
      set_face_level_next (jf3, level_face_new);

      ItChild it_child (rank,jf3);

      int jc3[3];
      while (it_child.next(jc3)) {

	Index index_child = index_.index_child(jc3);
	ItFace it_face_child = this->it_face(min_face_rank,index_child,jc3,jf3);

	int kf3[3];
	while (it_face_child.next()) {
	  it_face_child.face(kf3);
	  set_child_face_level_next(jc3,kf3,level_face_new);
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void Block::adapt_coarsen_()
{
#ifdef DEBUG_ADAPT
  CkPrintf ("%s COARSEN\n",name().c_str());
  fflush(stdout);
#endif
  const int level = this->level();
  
  // send data to parent

  ASSERT2 ("adapt_coarsen_()",
	  "Leaf = %s must be true and level = %d must be greater than 0",
	   is_leaf()?"true":"false", level,
	   is_leaf() && level > 0);

  Monitor * monitor = simulation()->monitor();
  if (monitor->is_verbose()) {
    char buffer [80];
    int v3[3];
    index().values(v3);
    sprintf (buffer,"Block %s (%x %x %x) is coarsening",name().c_str(),
	     v3[0],v3[1],v3[2]);
    monitor->print("Adapt",buffer);
  }

  // Create FieldFace for coarsening field data to parent

  int ic3[3];
  index_.child(level,&ic3[0],&ic3[1],&ic3[2]);
  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  std::vector<int> field_list;
  FieldFace * field_face = create_face
    (if3, ic3, lg3, refresh_coarse, field_list);

  const Index index_parent = index_.index_parent();

  // Create data message object to send

  DataMsg * data_msg = new DataMsg;

  // @@@ should be true but ~FieldFace() crashed in Sedov
  data_msg -> set_field_face (field_face,false);
  data_msg -> set_field_data (data()->field_data(),false);
  data_msg -> set_particle_data (data()->particle_data(),false);

  // copy face levels
  int nf = face_level_curr_.size();
  int face_level_curr[nf] = {0};
  
  for (int i=0; i<nf; i++) face_level_curr[i] = face_level_curr_[i];

  MsgCoarsen * msg = new MsgCoarsen (nf,face_level_curr,ic3);
  msg->set_data_msg (data_msg);

  thisProxy[index_parent].p_adapt_recv_child (msg);
  
}

//----------------------------------------------------------------------


void Block::p_adapt_recv_child (MsgCoarsen * msg)
{

  performance_start_(perf_adapt_update);

#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG p_refresh_store()\n",CkMyPe());
  fflush(stdout);
  if (msg->field_face()) msg->field_face()->print("called store");
#endif

  msg->update(data());

  int * ic3 = msg->ic3();
  int * child_face_level_curr = msg->face_level();

  // copy child face level array
  const int min_face_rank = 
    simulation()->config()->adapt_min_face_rank;
  Index index_child = index_.index_child(ic3);
  ItFace it_face_child = this->it_face(min_face_rank,index_child);
  int of3[3];

  while (it_face_child.next()) {
    it_face_child.face(of3);
    int level_child = child_face_level_curr[IF3(of3)];
    set_child_face_level_curr(ic3,of3,level_child);
  }

  // copy face level array
  ItFace it_face = this->it_face(min_face_rank,index_);
  while (it_face.next()) {
    it_face.face(of3);
    int level_child = child_face_level_curr[IF3(of3)];
    int opf3[3];
    if (parent_face_(opf3,of3,ic3)) {
      set_face_level_curr(opf3,level_child);
    }
  }

  // I am a leaf on the wind
  is_leaf_=true;

#ifdef DEBUG_ADAPT
  CkPrintf ("%s p_adapt_recv_child is_leaf <- 1\n",name().c_str());
#endif

  // Can now safely notify child to delete itself 
  adapt_delete_child_(index_child);

  age_ = 0;

  delete msg;

  performance_stop_(perf_adapt_update);
  performance_start_(perf_adapt_update_sync);
}


//----------------------------------------------------------------------

//----------------------------------------------------------------------

void Block::p_adapt_delete()
{
  performance_start_(perf_adapt_end);
#ifdef DEBUG_ADAPT
  CkPrintf ("%s DELETING\n",name().c_str());
#endif
  delete_ = true;
  performance_stop_(perf_adapt_end);
  performance_start_(perf_adapt_end_sync);
}

//======================================================================

void Block::initialize_child_face_levels_()
{
  const int  rank         = this->rank();
  const int level = this->level();

  if (level < 0) return;

  int ic3[3];
  ItChild it_child(rank);

  const int min_face_rank = 
    simulation()->config()->adapt_min_face_rank;

  // For each child

  while (it_child.next(ic3)) {

    // For each child face

    Index index_child = index_.index_child(ic3);
    ItFace it_face = this->it_face(min_face_rank,index_child);
    int if3[3];
    while (it_face.next()) {
      it_face.face(if3);
      int ip3[3];
      parent_face_(ip3,if3,ic3);
      Index in = neighbor_ (if3,&index_child);
      Index inp = in.index_parent();
      // Determine level for the child's face
      int level_child_face = (inp == thisIndex) ? 
	level + 1 : face_level(ip3);
      set_child_face_level_curr(ic3,if3, level_child_face);
    }

    if3[0]=if3[1]=if3[2]=0;
    set_child_face_level_curr(ic3,if3,level+1);

  }

  child_face_level_next_ = child_face_level_curr_;
}

//----------------------------------------------------------------------

bool Block::parent_face_
(int       ip3[3],
 const int if3[3],
 const int ic3[3]) const
{ 
  bool retval = false;
  for (int i=0; i<3; i++) {
    const bool inner_face_lower = (if3[i] == +1 && ic3[i] == 0);
    const bool inner_face_upper = (if3[i] == -1 && ic3[i] == 1);
    ip3[i] = (inner_face_lower || inner_face_upper) ? 0 : if3[i];

    retval = (retval || if3[i]);
  }
  return retval;
}

