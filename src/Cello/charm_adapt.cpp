
// See LICENSE_CELLO file for license and copyright information

/// @file     charm_new_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions
///
///----------------------------------------------------------------------

const char * adapt_str[] = {"unknown","coarsen","same","refine"};

// #define DEBUG_ADAPT

//--------------------------------------------------

static char buffer [256];

#ifdef DEBUG_ADAPT

#define PUT_LEVEL(INDEX,IC3,IF3,LEVEL_NOW,LEVEL_NEW,MSG)	\
  {								\
    std::string bit_str = INDEX.bit_string(-1,2);		\
    sprintf (buffer,"%12s %8s "					\
	     "if3 %2d %2d %2d  ic3 %d %d %d  "			\
	     "%d -> %d :%d",					\
	     MSG,bit_str.c_str(),				\
	     IF3[0],IF3[1],IF3[2],				\
	     IC3[0],IC3[1],IC3[2],				\
	     LEVEL_NOW,LEVEL_NEW,__LINE__);			\
    index_.print(buffer,-1,2,false,simulation());		\
    check_child_(IC3,"PUT_LEVEL",__FILE__,__LINE__);		\
    check_face_(IF3,"PUT_LEVEL",__FILE__,__LINE__);		\
    thisProxy[INDEX].p_adapt_recv_neighbor_level				\
      (index_,IC3,IF3,LEVEL_NOW,LEVEL_NEW);			\
  }
#else /* DEBUG_ADAPT */

#define PUT_LEVEL(INDEX,IC3,IF3,LEVEL_NOW,LEVEL_NEW,MSG)		\
  thisProxy[INDEX].p_adapt_recv_neighbor_level (index_,IC3,IF3,LEVEL_NOW,LEVEL_NEW);
#endif /* DEBUG_ADAPT */
//--------------------------------------------------

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//======================================================================

void CommBlock::adapt_begin_()
{

  set_leaf();

  const int rank = simulation()->dimension();
  sync_coarsen_.set_stop(NC(rank));
  sync_coarsen_.clear();

  const int initial_cycle = simulation()->config()->initial_cycle;
  const bool is_first_cycle = (initial_cycle == cycle());

  int level_maximum = is_first_cycle ? 
    simulation()->config()->initial_max_level :
    simulation()->config()->mesh_max_level;

  level_new_ = adapt_compute_desired_level_(level_maximum);

  control_sync (phase_sync_adapt_called);
}

//----------------------------------------------------------------------

bool CommBlock::do_adapt_()
{

  int adapt_interval = simulation()->config()->mesh_adapt_interval;

  return ((adapt_interval && ((cycle_ % adapt_interval) == 0)));

}

//----------------------------------------------------------------------

void CommBlock::adapt_called_()

{
  if (is_leaf()) {
    adapt_send_neighbors_levels(level_new_);
  }

  control_sync (phase_sync_adapt_next);

}

//----------------------------------------------------------------------

int CommBlock::adapt_compute_desired_level_(int level_maximum)
{

  if (! is_leaf()) return adapt_same;

  adapt_ = adapt_unknown;

  int level = this->level();
  int level_desired = level;

  FieldDescr * field_descr = simulation()->field_descr();

  int index_refine = 0;

  Problem * problem = simulation()->problem();
  Refine * refine;

  while ((refine = problem->refine(index_refine++))) {

    adapt_ = std::max(adapt_,refine->apply(this, field_descr));

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

void CommBlock::adapt_next_()
{

#ifdef DEBUG_ADAPT
  if (level() != level_new_) {
    sprintf (buffer,"is leaf %d level %d -> %d",is_leaf(),level(),level_new_);
    index_.print(buffer,-1,2,false,simulation());
  }
#endif

  update_levels_();

  if (is_leaf()) {
#ifdef DEBUG_ADAPT
    // sprintf (buffer,"level %d level_new_ %d",level(),level_new_);
    // index_.print(buffer,-1,2,false,simulation());
#endif
    if (level() < level_new_) adapt_refine_();
    if (level() > level_new_) adapt_coarsen_();
  }

  control_sync (phase_sync_adapt_end);
}

//----------------------------------------------------------------------

void CommBlock::adapt_end_()
{
  if (delete_) {

    // --------------------------------------------------
    // ENTRY: #6 SimulationCharm::adapt_exit_() -> ckDestroy()
    // ENTRY: if delete
    // ENTRY: adapt phase
    // --------------------------------------------------
    ckDestroy();
    // --------------------------------------------------

    return;
  }

  //  if (index_.is_root()) {
  //    thisProxy.doneInserting();
  //  }

  set_leaf();

  adapt_exit_();
}

//----------------------------------------------------------------------

void CommBlock::adapt_refine_()
{

  adapt_ = adapt_unknown;

  const int rank = simulation()->dimension();
  
  int nx,ny,nz;
  block()->field_block()->size(&nx,&ny,&nz);

  // For each new child

  int ic3[3];
  ItChild it_child (rank);
  while (it_child.next(ic3)) {

    Index index_child = index_.index_child(ic3);

    // If child doesn't exist yet

    if ( ! is_child_(index_child) ) {

      // Prolong data

      int narray = 0;  
      char * array = 0;
      int iface[3] = {0,0,0};
      bool lghost[3] = {true,true,true};
      FieldFace * field_face = 
	load_face_ (&narray,&array, iface,ic3,lghost, op_array_prolong);

      int num_field_blocks = 1;
      bool testing = false;

      const Factory * factory = simulation()->factory();

      // Create child block

      factory->create_block 
	(&thisProxy, index_child,
	 nx,ny,nz,
	 num_field_blocks,
	 adapt_step_ + 1,
	 cycle_,time_,dt_,
	 narray, array, op_array_prolong,
	 27,&child_face_level_[27*IC3(ic3)],
	 testing,
	 simulation());

      delete field_face;

      children_.push_back(index_child);

    }
  }
  set_leaf();
}

//----------------------------------------------------------------------

void CommBlock::adapt_delete_child_(Index index_child)
{
  // --------------------------------------------------
  // ENTRY: #9 CommBlock::adapt_delete_child_()-> CommBlock::p_adapt_delete()
  // ENTRY: child block
  // ENTRY: adapt phase
  // --------------------------------------------------
  thisProxy[index_child].p_adapt_delete();
  // --------------------------------------------------

  if (sync_coarsen_.next()) {
    children_.clear();
  }
}

//----------------------------------------------------------------------

void CommBlock::adapt_send_neighbors_levels(int level_new)
{
  const int level        = this->level();
  const int rank         = simulation()->dimension();
  const int rank_refresh = simulation()->config()->field_refresh_rank;

  ItFace it_face(rank,rank_refresh);
  int of3[3];

  while (it_face.next(of3)) {

    int ic3[3] = {0,0,0};

    Index index_neighbor = neighbor_(of3);

    const int level_face = face_level (of3);

    if (level_face == level) {

      // SEND-SAME: Face and level are sent to unique
      // neighboring block in the same level

      // --------------------------------------------------
      // ENTRY: #10 CommBlock::adapt_send_neighbors_levels()-> CommBlock::p_adapt_recv_neighbor_level()
      // ENTRY: same-level neighbor
      // ENTRY: adapt phase
      // --------------------------------------------------
      thisProxy[index_neighbor].p_adapt_recv_neighbor_level 
	(index_,ic3,of3,level,level_new);
      //      PUT_LEVEL (index_neighbor,ic3,of3,level,level_new,"send-same");
      // --------------------------------------------------

    } else if (level_face == level - 1) {

      // SEND-COARSE: Face, level, and child indices are sent to
      // unique neighboring block in the next-coarser level

      index_.child (level,&ic3[0],&ic3[1],&ic3[2]);

      int op3[3];
      parent_face_(op3,of3,ic3);

      // avoid redundant calls
      if (op3[0]==of3[0] && 
	  op3[1]==of3[1] && 
	  op3[2]==of3[2]) {

	Index index_uncle = index_neighbor.index_parent();

	// --------------------------------------------------
	// ENTRY: #11 CommBlock::adapt_send_neighbors_levels()-> CommBlock::p_adapt_recv_neighbor_level()
	// ENTRY: coarse-level neighbor
	// ENTRY: adapt phase
	// --------------------------------------------------
	PUT_LEVEL (index_uncle,ic3,of3,level,level_new,"send-coarse");
	// --------------------------------------------------

      }

    } else if (level_face == level + 1) {

      // SEND-FINE: Face and level are sent to all nibling
      // blocks in the next-finer level along the face.

      const int if3[3] = {-of3[0],-of3[1],-of3[2]};
      ItChild it_child(rank,if3);
      while (it_child.next(ic3)) {
	Index index_nibling = index_neighbor.index_child(ic3);

	// --------------------------------------------------
	// ENTRY: #12 CommBlock::adapt_send_neighbors_levels()-> CommBlock::p_adapt_recv_neighbor_level()
	// ENTRY: fine-level neighbor
	// ENTRY: adapt phase
	// --------------------------------------------------
	PUT_LEVEL (index_nibling,ic3,of3,level,level_new,"send-fine");
	// --------------------------------------------------

      }

    } else {
      std::string bit_str = index_.bit_string(-1,2);
      WARNING3 ("CommBlock::notify_neighbor()",
		"%s level %d and face level %d differ by more than 1",
		bit_str.c_str(),level,level_face);
    }

  }
}

//----------------------------------------------------------------------

void CommBlock::p_adapt_recv_neighbor_level
(
 Index index_caller,   // index of the calling neighbor
 int ic3[3],
 int if3[3], 
 int level_face,
 int level_face_new
 )
{
  int level_new = level_new_;

  const int level        = this->level();

  const int of3[3] = {-if3[0],-if3[1],-if3[2] };
  const int rank = simulation()->dimension();
  const int rank_refresh = simulation()->config()->field_refresh_rank;

  if (is_leaf()) {

    if (level == level_face) {

      // RECV-SAME: Face and level are received from unique
      // neighbor.  Unique face level is updated, and levels on
      // possibly multiple faces of multiple children are updated.

      set_face_level_new (of3, level_face_new);

      ItChild it_child (rank,of3);
      int jc3[3];
      while (it_child.next(jc3)) {

	int jf3[3];

	ItFace it_face (rank,rank_refresh,jc3,of3);

	while (it_face.next(jf3)) {
	  set_child_face_level_new(jc3,jf3,level_face_new);
	}
      }

    } else if (level == level_face - 1) {

      // RECV-FINE: Face, level, and sender child indices are received
      // from possibly non-unique neighbor for the corresponding face.
      // If neighbor level indicates the neighbor needs to refine,
      // then refinement is triggered in this block as well.  Face
      // levels are ignored, since the face either stays the same, or
      // refinement is triggered, in which case the child face levels
      // are used.  The neighbor level for the unique face and unique
      // child facing the sending child is updated.

      set_face_level_new (of3, level_face_new);

      ItChild it_child (rank,of3);
      int jc3[3];
      Index index_neighbor = neighbor_(of3).index_child(ic3);
      while (it_child.next(jc3)) {

	Index index_child = index_.index_child(jc3);

	int jf3[3];
	ItFace it_face (rank,rank_refresh,jc3,of3);

	while (it_face.next(jf3)) {

	  Index in = neighbor_(jf3,&index_child);

	  if (in == index_neighbor) {
	    set_child_face_level_new(jc3,jf3,level_face_new);
	  }
	}
       
      }

    } else if (level == level_face + 1) {

      // RECV-COARSE: Face and level are received from unique
      // neighbor.  Possibly multiple faces of block are updated
      // corresponding to the coarse neighbor's face.  Levels of
      // possibly multiple faces of possibly multiple child faces are
      // updated.

      ItFace it_face (rank,rank_refresh,ic3,of3);

      int jf3[3];
      while (it_face.next(jf3)) {
	set_face_level_new (jf3, level_face_new);

	ItChild it_child (rank,jf3);

	int jc3[3];
	while (it_child.next(jc3)) {

	  int kf3[3];

	  ItFace it_face (rank,rank_refresh,jc3,jf3);

	  while (it_face.next(kf3)) {
	    set_child_face_level_new(jc3,kf3,level_face_new);
	  }
	}
      }

    } else  {

      WARNING2 ("CommBlock::notify_neighbor()",
		"level %d and face level %d differ by more than 1",
		level,level_face);
    }

    // Don't coarsen if any siblings don't coarsen

      bool is_sibling = (index_.level() > 0 && index_caller.level() > 0) ? (index_caller.index_parent() == index_.index_parent()) : false;

    bool is_nephew = (index_.level() > 0 && index_caller.level() > 1) ? (index_caller.index_parent().index_parent() == index_.index_parent()) : false;
    bool is_coarsening = level_new < level;
    bool is_finer_neighbor = level_face_new > level_new;

    if (is_coarsening && ((is_sibling && is_finer_neighbor) || is_nephew )) {

#ifdef DEBUG_ADAPT
      // index_.print("DEBUG not coarsening",-1,2,false,simulation());
#endif

      level_new = level;
    }

    if (level_new < level_face_new - 1) {

      // restrict new level to within 1 of neighbor
      level_new = level_face_new - 1;

    }

    // notify neighbors if level_new has changed

    if (level_new != level_new_) {
      ASSERT2 ("CommBlock::p_adapt_recv_neighbor_level()",
	       "level_new %d level_new_ %d\n", level_new,level_new_,
	       level_new > level_new_);
      level_new_ = level_new;
      adapt_send_neighbors_levels(level_new_);
    }

  } else { // not a leaf

    Index index_neighbor = neighbor_(of3);

    // Forward to children if internal node
    
    ItChild it_child (rank,if3);
    int jc3[3];
    while (it_child.next(jc3)) {
      Index index_nibling = index_neighbor.index_child(jc3);
      
      //      index_.print("p_adapt_recv_neighbor_level() recurse",-1,2,false,simulation());
      // --------------------------------------------------
      // ENTRY: #13 CommBlock::p_adapt_recv_neighbor_level()-> CommBlock::p_get_neighbor_level()
      // ENTRY: non-leaf recurse on children
      // ENTRY: adapt phase
      // --------------------------------------------------
      PUT_LEVEL (index_nibling,ic3,if3,level_face,level_face_new,"RECURSE");
      // --------------------------------------------------

    }
  }
}

//----------------------------------------------------------------------

void CommBlock::adapt_coarsen_()
{

  const int level = this->level();
  
  // send data to parent

  if (level > 0 && is_leaf()) {

    Index index_parent = index_.index_parent();
    int ic3[3] = {1,1,1};
    index_.child(level,&ic3[0],&ic3[1],&ic3[2]);

    // copy block data
    int narray; 
    char * array;
    int iface[3] = {0,0,0};
    bool lghost[3] = {true,true,true};
    FieldFace * field_face = 
      load_face_(&narray,&array, iface, ic3, lghost, op_array_restrict);

    // copy face levels
    int nf = face_level_.size();
    int face_level[nf];
    for (int i=0; i<nf; i++) face_level[i] = face_level_[i];

    // send child data to parent

    // --------------------------------------------------
    // ENTRY: #14 CommBlock::adapt_coarsen_()-> CommBlock::p_adapt_send_child_data()
    // ENTRY: parent if leaf() and level > 0
    // ENTRY: adapt phase
    // --------------------------------------------------
    thisProxy[index_parent].p_adapt_send_child_data
      (ic3, narray,array, nf,face_level);
    // --------------------------------------------------

    delete field_face;

  }
}

//----------------------------------------------------------------------

void CommBlock::p_adapt_send_child_data
(
 int ic3[3],
 int na, char * array,
 int nf, int * child_face_level
 )
{
  // copy array
  int iface[3] = {0,0,0};
  bool lghost[3] = {true,true,true};
  store_face_(na,array, iface, ic3, lghost, op_array_restrict);

  // copy child face level and face level
  const int rank = simulation()->dimension();
  const int  rank_refresh = simulation()->config()->field_refresh_rank;
  ItFace  it_face(rank,rank_refresh);
  int of3[3];
  while (it_face.next(of3)) {
    int level_child = child_face_level[IF3(of3)];
    set_child_face_level(ic3,of3,level_child);
    int opf3[3];
    if (parent_face_(opf3,of3,ic3)) {
      set_face_level(opf3,level_child);
    }
  }

  Index index_child = index_.index_child(ic3);

  adapt_delete_child_(index_child);

}

//----------------------------------------------------------------------

void CommBlock::p_adapt_delete()
{
  delete_ = true;
    // // --------------------------------------------------
    // // ENTRY: #6 SimulationCharm::adapt_exit_() -> ckDestroy()
    // // ENTRY: if delete
    // // ENTRY: adapt phase
    // // --------------------------------------------------
    // ckDestroy();
    // // --------------------------------------------------
}

//======================================================================

void CommBlock::initialize_child_face_levels_()
{
  const int  rank         = simulation()->dimension();
  const int  rank_refresh = simulation()->config()->field_refresh_rank;
  const int level = this->level();

  int ic3[3];
  ItChild it_child(simulation()->dimension());

  // For each child

  while (it_child.next(ic3)) {

    Index index_child = index_.index_child(ic3);
    int if3[3];
    ItFace it_face(rank,rank_refresh);

    // For each child face

    while (it_face.next(if3)) {

      int ip3[3];
      parent_face_(ip3,if3,ic3);
      Index in = neighbor_ (if3,&index_child);
      Index inp = in.index_parent();
      // Determine level for the child's face
      int level_child_face = (inp == thisIndex) ? 
	level + 1 : face_level(ip3);
      // Set the child's face to the level
      set_child_face_level(ic3,if3, level_child_face);
    }

    if3[0]=if3[1]=if3[2]=0;
    set_child_face_level(ic3,if3,level+1);

  }

  child_face_level_new_ = child_face_level_;
}

//----------------------------------------------------------------------

bool CommBlock::parent_face_
(int       ip3[3],
 const int if3[3],
 const int ic3[3]) const
{
  ip3[0] = if3[0];
  ip3[1] = if3[1];
  ip3[2] = if3[2];
  if (if3[0] == +1 && ic3[0] == 0) ip3[0] = 0;
  if (if3[1] == +1 && ic3[1] == 0) ip3[1] = 0;
  if (if3[2] == +1 && ic3[2] == 0) ip3[2] = 0;
  if (if3[0] == -1 && ic3[0] == 1) ip3[0] = 0;
  if (if3[1] == -1 && ic3[1] == 1) ip3[1] = 0;
  if (if3[2] == -1 && ic3[2] == 1) ip3[2] = 0;
  return (if3[0] || if3[1] || if3[2]);
}

//----------------------------------------------------------------------

FieldFace * CommBlock::load_face_
(
 int *   n, char ** a,
 int if3[3], int ic3[3], bool lg3[3],
 int op_array_type
 )
{
  FieldFace * field_face = create_face_ (if3,ic3,lg3, op_array_type);

  field_face->load(n, a);
  return field_face;
}

//----------------------------------------------------------------------

void CommBlock::store_face_
(
 int n, char * a, 
 int if3[3], int ic3[3], bool lg3[3],
 int op_array_type
 )
{
  FieldFace * field_face = create_face_ (if3,ic3,lg3, op_array_type);

  field_face->store(n, a);
  delete field_face;
}

//----------------------------------------------------------------------

FieldFace * CommBlock::create_face_
(int if3[3], int ic3[3], bool lg3[3],
 int op_array_type
 )
{
  Problem * problem        = simulation()->problem();
  FieldDescr * field_descr = simulation()->field_descr();
  FieldBlock * field_block = block_->field_block();

  FieldFace * field_face = new FieldFace (field_block,field_descr);

  if (op_array_type == op_array_restrict) {

    field_face->set_restrict(problem->restrict(),ic3[0],ic3[1],ic3[2]);

  } else if (op_array_type == op_array_prolong) {

    field_face->set_prolong(problem->prolong(),  ic3[0],ic3[1],ic3[2]);

  }

  field_face->set_face (if3[0],if3[1],if3[2]);
  field_face->set_ghost(lg3[0],lg3[1],lg3[2]);

  return field_face;
}

