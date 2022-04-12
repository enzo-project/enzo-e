// See LICENSE_CELLO file for license and copyright information

/// @file     control_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions.
/// @ingroup  Control
///
/// This file controls adaptive mesh refinement on a distributed
/// array of octrees.

//--------------------------------------------------

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define TRACE_ADAPT
// #define TRACE_REFINE
// #define TRACE_COARSEN
// #define TRACE_DONE_INSERTING
// #define DEBUG_ADAPT

// #define WRITE_NEIGHBORS

#define DEBUG_CYCLE_START 0

#ifdef TRACE_ADAPT
#   undef TRACE_ADAPT
#   define TRACE_ADAPT(MSG,BLOCK)                                       \
  if (BLOCK->cycle() >= DEBUG_CYCLE_START) {                            \
    CkPrintf ("TRACE_ADAPT %s %s step %d ready %d self conv %d neighb conv %d\n", \
              BLOCK->name().c_str(),                                    \
              std::string(MSG).c_str(),                                 \
              BLOCK->adapt_step_,                                       \
              BLOCK->adapt_ready_?1:0 ,                                 \
              BLOCK->adapt_.is_converged()?1:0,                         \
              BLOCK->adapt_.neighbors_converged()?1:0);                 \
    fflush(stdout);                                                     \
  }
#else
#   define TRACE_ADAPT(MSG,BLOCK) /* ... */
#endif

//======================================================================

/// @brief First function in the adapt phase: apply local refinement criteria.
///
/// adapt_begin_() computes the desired refinement level for the block
/// using adapt_compute_desired_level_(), after which it calls
/// adapt_called_() with nearest-neighbor synchronization.

void Block::adapt_enter_()
{
  TRACE_ADAPT("adapt_enter_",this);
  if ( do_adapt_()) {

    adapt_begin_();

  } else {

    adapt_exit_();

  }
}

//----------------------------------------------------------------------

void Block::adapt_begin_()
{
  TRACE_ADAPT("adapt_begin_",this);

  cello::simulation()->set_phase(phase_adapt);

  if (is_leaf()) {

#ifdef WRITE_NEIGHBORS
    adapt_.write("adapt",this,DEBUG_CYCLE_START);
#endif

  // Evaluate local mesh refinement criteria
    const int level_maximum = cello::config()->mesh_max_level;
    level_next_ = adapt_compute_desired_level_(level_maximum);

    // Reset adapt level bounds for next adapt phase
    adapt_.reset_bounds();
    adapt_.initialize_self(index_,level_next_,index_.level());
    adapt_.update_bounds();
  }
#ifdef DEBUG_ADAPT
  CkPrintf ("DEBUG_ADAPT %s level_next = %d\n",name().c_str(),level_next_);
#endif
  const int min_face_rank = cello::config()->adapt_min_face_rank;
  control_sync_neighbor (CkIndex_Block::p_adapt_called(),
			 sync_id_adapt_begin,
			 min_face_rank,
			 neighbor_leaf,0);
}

//----------------------------------------------------------------------

/// @brief Second step of the adapt phase: tell neighbors desired level.
///
/// Call adapt_send_level() to send neighbors desired
/// levels, after which adapt_next_() is called with quiescence
/// detection.
void Block::adapt_called_()
{
  TRACE_ADAPT("adapt_called_",this);
  if (! is_leaf()) {
    TRACE_ADAPT("adapt_barrier [not leaf]",this);
    adapt_barrier_();
    return;
  }
  adapt_send_level();

}

//----------------------------------------------------------------------
void Block::adapt_barrier_()
{
  if (! adapt_balanced_) {
    adapt_balanced_ = true;
    TRACE_ADAPT("calling contribute",this);
    int changed = (is_leaf() && (level() != level_next_)) ? 1 : 0;
    CkCallback callback = CkCallback
      (CkIndex_Block::r_adapt_next(nullptr), 
       proxy_array());
    adapt_ready_ = true;
    contribute(sizeof(int),&changed,CkReduction::sum_int, callback);
  }
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
  TRACE_ADAPT("adapt_next contribute called_",this);
  update_levels_();

  level_next_ = adapt_.level_min();
  
  if (is_leaf()) {
    if (level() < level_next_) {
#ifdef TRACE_REFINE
      CkPrintf ("TRACE_REFINE %s\n",name().c_str());
#endif
      adapt_refine_();
    }
    if (level() > level_next_) {
#ifdef TRACE_COARSEN
      CkPrintf ("TRACE_COARSEN %s\n",name().c_str());
#endif
      adapt_coarsen_();
    }
  }
  control_sync_quiescence (CkIndex_Main::p_adapt_update());

}

//----------------------------------------------------------------------
void Block::update_levels_ ()
{
  if (! is_leaf()) return;
  adapt_.update_curr_from_next();
  child_face_level_curr_ = child_face_level_next_;

  std::vector<Index> index_list;
  // Save list of indices before updating them
  const int n = adapt_.num_neighbors();
  index_list.resize(n);
  for (int i=0; i<n; i++) {
    index_list[i] = adapt_.index(i);
  }
  // Loop through neighbor levels and refine/coarsen as needed
  for (int i=0; i<n; i++) {
    int level_min;
    int level_max;
    bool can_coarsen;
    adapt_.get_neighbor_level_bounds
      (index_list[i],&level_min,&level_max,&can_coarsen);
    if (level_min > index_list[i].level()) {
      // level_min is larger, must refine
      adapt_.refine_neighbor(index_list[i]);
    } else if (level_min < index_list[i].level() && can_coarsen) {
      // level_min is smaller and can coarsen
      adapt_.coarsen_neighbor(index_list[i]);
    }
  }
}

//----------------------------------------------------------------------

void Block::adapt_update_()
{
  TRACE_ADAPT("adapt_update_",this);
  if (index_.is_root() && (adapt_changed_ != 0)) {
#ifdef TRACE_DONE_INSERTING
    CkPrintf ("TRACE_DONE_INSERTING\n");
#endif
    thisProxy.doneInserting();
  }
  adapt_end_();
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
  TRACE_ADAPT("adapt_end_",this);
  adapt_.reset_face_level(Adapt::LevelType::last);

  sync_coarsen_.reset();
  sync_coarsen_.set_stop(cello::num_children());

  const int initial_cycle = cello::config()->initial_cycle;
  const bool is_first_cycle = (initial_cycle == cycle());
  const int level_maximum = cello::config()->mesh_max_level;

  bool adapt_again = (is_first_cycle && (adapt_step_ < level_maximum));
  adapt_step_++;
  adapt_ready_ = false;
  adapt_balanced_ = false;

  if (adapt_again) {
    control_sync_quiescence (CkIndex_Main::p_adapt_enter());
  } else {
    control_sync_quiescence (CkIndex_Main::p_adapt_exit());
  }

}

//----------------------------------------------------------------------

/// @brief Return whether the adapt phase should be called this cycle.
bool Block::do_adapt_()
{
  int adapt_interval = cello::config()->adapt_interval;

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
  TRACE_ADAPT("adapt_compute_desired_level_",this);
  if (! is_leaf()) return adapt_same;

  int adapt = adapt_unknown;

  int level = this->level();
  int level_desired = level;

  Problem * problem = cello::problem();
  Refine * refine;

  int index_refine = 0;
  while ((refine = problem->refine(index_refine++))) {

    Schedule * schedule = refine->schedule();

    if ((schedule==NULL) || schedule->write_this_cycle(cycle(),time()) ) {
      adapt = std::max(adapt,refine->apply(this));
    }

  }
  const int initial_cycle = cello::config()->initial_cycle;
  const bool is_first_cycle = (initial_cycle == cycle());

  if (adapt == adapt_coarsen && level > 0 && ! is_first_cycle)
    level_desired = level - 1;
  else if (adapt == adapt_refine  && level < level_maximum)
    level_desired = level + 1;
  else {
    adapt = adapt_same;
    level_desired = level;
  }

  return level_desired;
}

//----------------------------------------------------------------------

void Block::adapt_refine_()
{
  TRACE_ADAPT("adapt_refine",this);
  Monitor * monitor = cello::monitor();
  if (monitor->is_verbose()) {
    char buffer [80];
    int v3[3];
    index().values(v3);
    sprintf (buffer,"Block %s (%x %x %x) is refining",name().c_str(),
	     v3[0],v3[1],v3[2]);
    monitor->print("Adapt",buffer);
  }

  int nx,ny,nz;
  data()->field_data()->size(&nx,&ny,&nz);

  // First scatter particles to children to avoid multiple passes

  ParticleData * particle_list[8] = {0};

  ParticleDescr * p_descr = cello::simulation()->particle_descr();

  const int nc = cello::num_children();

  for (int i=0; i<nc; i++) {
    particle_list[i] = new ParticleData;
    particle_list[i]->allocate(p_descr);
  }

  Particle particle = data()->particle();

  // partition particles into 2x2x2 array to send to children

  particle_scatter_children_ (particle_list,particle);

  // For each new child

  const int rank = cello::rank();

  ItChild it_child (rank);
  int ic3[3];
  while (it_child.next(ic3)) {

    Index index_child = index_.index_child(ic3);

    // If child doesn't exist yet

    if ( ! is_child_(index_child) ) {

      // Create FieldFace for interpolating field data to child ic3[]

      int if3[3] = {0,0,0};
      int g3[3];
      cello::field_descr()->ghost_depth(0,g3,g3+1,g3+2);
      Refresh * refresh = new Refresh;
      refresh->add_all_data();
      FieldFace * field_face = create_face
	(if3,ic3,g3, refresh_fine, refresh, true);

      // Create data message object to send
      DataMsg * data_msg = new DataMsg;

      // @@@ should be true but ~FieldFace() crashes
      data_msg -> set_field_face (field_face,false);
      data_msg -> set_field_data (data()->field_data(),false);
      ParticleData * p_data = new ParticleData(*particle_list[IC3(ic3)]);
      data_msg -> set_particle_data (p_data,true);

      const Factory * factory = cello::simulation()->factory();

      // Create the child object with interpolated data

      int narray = 0;
      char * array = 0;
      int num_field_data = 1;

      factory->create_block
	(
	 data_msg,
	 thisProxy, index_child,
	 nx,ny,nz,
	 num_field_data,
	 adapt_step_,
	 cycle_,time_,dt_,
	 narray, array, refresh_fine,
	 27,
         &child_face_level_curr_.data()[27*IC3(ic3)],
         &adapt_,
	 cello::simulation());

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
  cello::simulation()->data_delete_particles(count);

  adapt_.set_valid(false);
  is_leaf_ = false;
  TRACE_ADAPT("adapt_refine exit",this);
}

//----------------------------------------------------------------------

void Block::particle_scatter_children_ (ParticleData * particle_list[],
					Particle particle)
{
  const int npa = cello::num_children();

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
      const bool * mask = nullptr;

      // ...extract particle position arrays
      std::vector<double> xa(np);
      std::vector<double> ya(np);
      std::vector<double> za(np);

      particle.position(it,ib,xa.data(),ya.data(),za.data());

      // ...and corresponding particle indices
      std::vector<int> index(np);

      if (is_float) {

	// absolute coordinates

	for (int ip=0; ip<np; ip++) {

	  const double x = xa[ip*d];
	  const double y = ya[ip*d];
	  const double z = za[ip*d];

	  const int rank = cello::rank();
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
      particle.scatter (it,ib, np, mask, index.data(), npa, particle_list);
      // ... delete scattered particles
      count += particle.delete_particles (it,ib,mask);

      delete [] mask;
    }
  }
  cello::simulation()->data_delete_particles(count);
}

//----------------------------------------------------------------------

void Block::adapt_delete_child_(Index index_child)
{
  thisProxy[index_child].p_adapt_delete();

  if (sync_coarsen_.next()) {
    children_.clear();
  }
}

//----------------------------------------------------------------------

void Block::adapt_send_level()
{

  adapt_ready_ = true;
  TRACE_ADAPT("adapt_send_level",this);
  if (!is_leaf()) return;
  const int level = this->level();

  int level_min;
  int level_max;
  bool can_coarsen;
  adapt_.update_bounds();
  adapt_.get_level_bounds(&level_min,&level_max,&can_coarsen);

  ItNeighbor it_neighbor = this->it_neighbor(index_);
  int of3[3];
  std::map<Index,int> index_count;
  std::map<Index,bool> index_first;
  while (it_neighbor.next(of3)) {
    Index index_neighbor = it_neighbor.index();
    index_first[index_neighbor] = true;
    ++index_count[index_neighbor];
  }

  std::map<Index,MsgAdapt *> msg_map;

  while (it_neighbor.next(of3)) {

    Index index_neighbor = it_neighbor.index();

    // Skip self if own neighbor (e.g. single-block periodic b.c.)
    if (index_neighbor == index_) continue;

    int ic3[3];
    it_neighbor.child(ic3);
    // int level_neighbor = index_neighbor.level();
    // if (level > level_neighbor) {
    //   index_.child(level,ic3,ic3+1,ic3+2);
    // } else if (level < level_neighbor) {
    //   index_neighbor.child(level_neighbor,ic3,ic3+1,ic3+2);
    // } else {
    //   ic3[0]=ic3[1]=ic3[2]=0;
    // }
    if (index_first[index_neighbor]) {

      msg_map[index_neighbor] = new MsgAdapt
        (adapt_step_,index_,ic3,of3,level,level_min,level_max,can_coarsen);
      index_first[index_neighbor] = false;
    } else {
      msg_map[index_neighbor]->add_face(of3);
    }

#ifdef DEBUG_ADAPT
    CkPrintf ("DEBUG_ADAPT %s : %s send_level of3 %d %d %dic3 %d %d %d\n",
              name().c_str(),name(index_neighbor).c_str(),
              of3[0],of3[1],of3[2],ic3[0],ic3[1],ic3[2]);
#endif
    if (--index_count[index_neighbor] == 0) {
#ifdef DEBUG_ADAPT
      CkPrintf ("DEBUG_ADAPT %s : %s send_level sending message\n",
                name().c_str(),name(index_neighbor).c_str());
#endif
      thisProxy[index_neighbor].p_adapt_recv_level (msg_map[index_neighbor]);
    }
  }
  TRACE_ADAPT("calling adapt_recv_level",this);
  adapt_recv_level();
}

void Block::p_adapt_recv_level (MsgAdapt * msg)
{
  if (!adapt_ready_) {
    // save message for later
    adapt_msg_list_.push_back(msg);
  } else {
    adapt_check_messages_();
    adapt_recv_level
      (
       msg->adapt_step_,
       msg->index_,
       msg->ic3_,
       msg->ofv_,
       msg->level_now_,
       msg->level_min_,
       msg->level_max_,
       msg->can_coarsen_);
    delete msg;
  }
}

void Block::adapt_recv_level()
{
  // Process any saved messages
  adapt_check_messages_();
  if (adapt_.neighbors_converged() && adapt_.is_converged()) {
    TRACE_ADAPT("adapt_barrier [self]",this);
    adapt_barrier_();
  }
}

void Block::adapt_check_messages_()
{
  for (int i=0; i<adapt_msg_list_.size(); i++) {
    MsgAdapt * msg = adapt_msg_list_[i];
    adapt_recv_level (
                      msg->adapt_step_,
                      msg->index_,
                      msg->ic3_,
                      msg->ofv_,
                      msg->level_now_,
                      msg->level_min_,
                      msg->level_max_,
                      msg->can_coarsen_);
  }
  adapt_msg_list_.clear();
}

void Block::adapt_recv_level
(
 int adapt_step,
 Index index_send,
 int ic3[3],
 std::vector<int> ofv[3],
 int level_face_curr,
 int level_face_new,
 int level_face_max,
 bool can_coarsen
 )
{
  bool changed = false;
  int level_min;
  for (int i=0; i<ofv[0].size(); i++) {

    int if3[3] = {ofv[0][i],ofv[1][i],ofv[2][i]};
#ifdef DEBUG_ADAPT
    CkPrintf ("DEBUG_ADAPT adapt_recv_level if3 %d %d %d\n",
              if3[0],if3[1],if3[2]);
#endif
    //  if (adapt_step != adapt_step_) return;

    ASSERT4("Block::p_adapt_recv_level()",
            "Mismatch in %s : %s adapt_recv_level actual %d != expected %d\n",
            name().c_str(),name(index_send).c_str(),adapt_step,adapt_step_,
            adapt_step == adapt_step_);
    TRACE_ADAPT("adapt_recv_level",this);
    int level_max;

    adapt_.update_neighbor
      (index_send,level_face_new,level_face_max,can_coarsen);
    changed = changed || adapt_.update_bounds();

    adapt_.get_level_bounds(&level_min,&level_max,&can_coarsen);

    ASSERT2("p_adapt_recv_level()",
            "Index level %d and indicated level %d mismatch",
            index_send.level(),level_face_curr,
            level_face_curr == index_send.level());
    performance_start_(perf_adapt_update);

    if (index_send.level() != level_face_curr) {
      PARALLEL_PRINTF
        ("%d level mismatch between index_send (%d) and level_face_curr (%d)",
         __LINE__,index_send.level(), level_face_curr);
      int nb3[3] = {2,2,2};
      index_.print("index_",-1,2,nb3,false,cello::simulation());
      index_send.print("index_",-1,2,nb3,false,cello::simulation());
    }

    adapt_.set_face_level_last
      (ic3,if3, level_face_new,level_face_max, can_coarsen);

    int level_next = level_next_;

    const int level = this->level();

    if ( ! is_leaf()) {

      ERROR1 ("Block::p_adapt_recv_level()",
              "Block %s is not a leaf",
              name().c_str());
    }

    const int of3[3] = {-if3[0],-if3[1],-if3[2] };

    ASSERT2("adapt_recv_level()",
            "level_face_curr %d and level %d differ by more than one",
            level_face_curr,level,
            std::abs(level_face_curr - level) <= 1);

    const int level_diff = (level_face_curr - level);

    adapt_recv(index_send,of3,ic3,level_face_new,level_diff,
               level_face_max,can_coarsen);

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

      // cannot if nephew not also coarsening
      //    if (is_nephew && (level_face_new > level)) {
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
  }
  if (changed) {
    level_next_ = level_min;
    adapt_send_level();
  }
  TRACE_ADAPT("testing convergence",this);
  if (adapt_.neighbors_converged() && adapt_.is_converged()) {
    TRACE_ADAPT("adapt_barrier [recv_level]",this);
    adapt_barrier_();
  }
  performance_stop_(perf_adapt_update);
  performance_start_(perf_adapt_update_sync);
}

//----------------------------------------------------------------------

void Block::adapt_recv
( Index index_send,
  const int of3[3],
  const int ic3[3], int level_face_new, int level_relative,
  int level_max, bool can_coarsen)
{

  const int rank = cello::rank();
  const int min_face_rank = cello::config()->adapt_min_face_rank;

  if (level_relative == 0 || level_relative == +1) {
    // RECV-SAME: Face and level are received from unique
    // neighbor.  Unique face level is updated, and levels on
    // possibly multiple faces of multiple children are updated.

    adapt_.set_face_level
      (of3, Adapt::LevelType::next,
       level_face_new,level_max,can_coarsen);

    ItChild it_child (rank,of3);
    int jc3[3];

    while (it_child.next(jc3)) {

      int kf3[3];

      Index index_child = index_.index_child(jc3);
      ItFace it_face = this->it_face(min_face_rank,index_child,jc3,of3);

      if (level_relative == 0) {
        while (it_face.next(kf3)) {
          set_child_face_level_next(jc3,kf3,level_face_new);
        }
      } else if (level_relative == +1) {
        while (it_face.next(kf3)) {

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
    while (it_face.next(jf3)) {

      adapt_.set_face_level (jf3, Adapt::LevelType::next,
                             level_face_new,level_max,can_coarsen);

      ItChild it_child (rank,jf3);

      int jc3[3];
      while (it_child.next(jc3)) {

        Index index_child = index_.index_child(jc3);
        ItFace it_face_child = this->it_face(min_face_rank,index_child,jc3,jf3);

        int kf3[3];
        while (it_face_child.next(kf3)) {
          set_child_face_level_next(jc3,kf3,level_face_new);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void Block::adapt_coarsen_()
{
  const int level = this->level();

  // send data to parent

  ASSERT2 ("adapt_coarsen_()",
	  "Leaf = %s must be true and level = %d must be greater than 0",
	   is_leaf()?"true":"false", level,
	   is_leaf() && level > 0);

  Monitor * monitor = cello::monitor();

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
  int g3[3] = {0,0,0};
  Refresh * refresh = new Refresh;
  refresh->add_all_data();

  FieldFace * field_face = create_face
    (if3, ic3, g3, refresh_coarse, refresh, true);

  const Index index_parent = index_.index_parent();

  // Create data message object to send

  DataMsg * data_msg = new DataMsg;

  // @@@ should be true but ~FieldFace() crashed in Sedov
  data_msg -> set_field_face (field_face,false);
  data_msg -> set_field_data (data()->field_data(),false);
  data_msg -> set_particle_data (data()->particle_data(),false);

  const int nf = adapt_.size_face_level(Adapt::LevelType::curr);

  MsgCoarsen * msg = new MsgCoarsen
    (nf,adapt_.vector_face_level(Adapt::LevelType::curr),ic3,&adapt_);

  msg->set_data_msg (data_msg);

  thisProxy[index_parent].p_adapt_recv_child (msg);

}

//----------------------------------------------------------------------


void Block::p_adapt_recv_child (MsgCoarsen * msg)
{
  TRACE_ADAPT("p_adapt_recv_child",this);

  performance_start_(perf_adapt_update);
  msg->update(data());
  int * ic3 = msg->ic3();
  int * child_face_level_curr = msg->face_level();
  // copy child face level array
  const int min_face_rank = cello::config()->adapt_min_face_rank;
  Index index_child = index_.index_child(ic3);
  ItFace it_face_child = this->it_face(min_face_rank,index_child);
  int of3[3];

  adapt_.coarsen(*msg->adapt_child_);

  while (it_face_child.next(of3)) {
    int level_child = child_face_level_curr[IF3(of3)];
    set_child_face_level_curr(ic3,of3,level_child);
  }

  // copy face level array
  ItFace it_face = this->it_face(min_face_rank,index_);
  while (it_face.next(of3)) {
    int level_child = child_face_level_curr[IF3(of3)];
    int opf3[3];
    if (parent_face_(opf3,of3,ic3)) {
      adapt_.set_face_level
        (opf3,Adapt::LevelType::curr,level_child);
    }
  }

  // I am a leaf on the wind
  is_leaf_=true;

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
  ckDestroy();
}

//======================================================================

void Block::initialize_child_face_levels_()
{
  const int  rank = cello::rank();
  const int level = this->level();

  if (level < 0) return;

  int ic3[3];
  ItChild it_child(rank);

  const int min_face_rank = cello::config()->adapt_min_face_rank;

  // For each child

  while (it_child.next(ic3)) {

    // For each child face

    Index index_child = index_.index_child(ic3);
    ItFace it_face = this->it_face(min_face_rank,index_child);
    int if3[3];
    while (it_face.next(if3)) {
      int ip3[3];
      parent_face_(ip3,if3,ic3);
      Index in = neighbor_ (if3,&index_child);
      Index inp = in.index_parent();
      // Determine level for the child's face
      const int level_face =
        adapt_.face_level(ip3,Adapt::LevelType::curr);
      int level_child_face = (inp == thisIndex) ?
        (level + 1) : level_face;
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
