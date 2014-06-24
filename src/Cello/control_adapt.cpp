// See LICENSE_CELLO file for license and copyright information

/// @file     control_adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions
/// @ingroup  Control

/// This file controls adaptive mesh refinement on a distributed
/// forest of octrees.
///
/// adapt_begin_()

//--------------------------------------------------
#define DEBUG_ADAPT
//--------------------------------------------------

#ifdef DEBUG_ADAPT

static char buffer [256];

#   define PUT_LEVEL(INDEX_SEND,INDEX_RECV,IC3,IF3,LEVEL_NOW,LEVEL_NEW,MSG)	\
  {								\
    std::string bit_str = INDEX_RECV.bit_string(-1,2);		\
    sprintf (buffer,"%12s %8s"					\
	     "if3 %2d %2d %2d  ic3 %d %d %d  "			\
	     "%d -> %d :%d",					\
	     MSG,bit_str.c_str(),				\
	     IF3[0],IF3[1],IF3[2],				\
	     IC3[0],IC3[1],IC3[2],				\
	     LEVEL_NOW,LEVEL_NEW,__LINE__);			\
    INDEX_SEND.print(buffer,-1,2,false,simulation());		\
    check_child_(IC3,"PUT_LEVEL",__FILE__,__LINE__);		\
    check_face_(IF3,"PUT_LEVEL",__FILE__,__LINE__);		\
    thisProxy[INDEX_RECV].p_adapt_recv_level				\
      (INDEX_SEND,IC3,IF3,LEVEL_NOW,LEVEL_NEW);			\
  }
#else /* DEBUG_ADAPT */

#   define PUT_LEVEL(INDEX_SEND,INDEX_RECV,IC3,IF3,LEVEL_NOW,LEVEL_NEW,MSG)	\
  thisProxy[INDEX_RECV].p_adapt_recv_level (INDEX_SEND,IC3,IF3,LEVEL_NOW,LEVEL_NEW);
#endif /* DEBUG_ADAPT */


#ifdef DEBUG_ADAPT

#ifdef CELLO_TRACE
#   define trace(A) fprintf (simulation()->fp_debug(),"%s:%d %s TRACE %s\n",__FILE__,__LINE__,name_.c_str(),A)
#else
#   define trace(A) /*  NULL */
#endif

#else

#   define trace(A) /*  NULL */

#endif
//--------------------------------------------------

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//======================================================================

/// @brief First function in the adapt phase: apply local refinement criteria.
/// 
/// adapt_begin_() computes the local desired refinement level using
/// adapt_compute_desired_level_(), after which it calls
/// adapt_called_() with nearest-neighbor synchronization.

void CommBlock::adapt_begin_()

{
  trace("adapt_begin 1");

  int level_maximum = simulation()->config()->mesh_max_level;

  level_new_ = adapt_compute_desired_level_(level_maximum);

  control_sync (sync_adapt_called,"neighbor");
}

//----------------------------------------------------------------------

/// @brief Second step of the adapt phase: tell neighbors desired level.
///
/// Call adapt_send_neighbors_levels() to send neighbors desired
/// levels, after which adapt_next_() is called with quiescence
/// detection.
void CommBlock::adapt_called_()
{
  trace("adapt_called 2");

  adapt_send_neighbors_levels(level_new_);

  control_sync (sync_adapt_next,"quiescence");

}

//----------------------------------------------------------------------

/// @brief Third step of the adapt phase: coarsen or refine according
/// to desired level.
///
/// Call update_levels_() to finalize face and child face levels,
/// then, if a leaf, refine or coarsen according to desired level
/// determined in adapt_called_().  Afterward, all CommBlocks call
/// adapt_end_().
void CommBlock::adapt_next_()
{

  trace("adapt_next 3");

#ifdef DEBUG_ADAPT
  if (level() != level_new_) {
    sprintf (buffer,"is leaf %d level %d -> %d",is_leaf(),level(),level_new_);
    index_.print(buffer,-1,2,false,simulation());
  }
#endif

  update_levels_();

  if (is_leaf()) {
    if (level() < level_new_) adapt_refine_();
    if (level() > level_new_) adapt_coarsen_();
  }

  control_sync (sync_adapt_end,"quiescence");
}

//----------------------------------------------------------------------

/// @brief Fourth step of the adapt phase: delete self if CommBlock has
/// been coarsened
///
/// This step deletes itself if it has been coarsened in this adapt
/// phase, then exits the adapt phase by directly calling adapt_exit_().
/// This is a separate phase since the quiescence call of this function
/// from the previous adapt_next_() step includes CommBlock's that have
/// been deleted.  
void CommBlock::adapt_end_()
{
  trace("adapt_next 4");

  if (index_.is_root()) {
    thisProxy.doneInserting();
  }

  if (delete_) {

#ifdef DEBUG_ADAPT
    index_.print("DESTROY",-1,2,false,simulation());
#endif
    ckDestroy();

    return;
  } else {

    for (size_t i=0; i<face_level_last_.size(); i++)
      face_level_last_[i] = 0;

    const int rank = simulation()->dimension();
    sync_coarsen_.set_stop(NC(rank));
    sync_coarsen_.clear();

    control_sync(sync_adapt_exit,"quiescence");
  }
}

//----------------------------------------------------------------------

/// @brief Return whether the adapt phase should be called this cycle.
bool CommBlock::do_adapt_()
{

  int adapt_interval = simulation()->config()->mesh_adapt_interval;

  return ((adapt_interval && ((cycle_ % adapt_interval) == 0)));

}

//----------------------------------------------------------------------

/// @brief Determine whether this CommBlock should refine, coarsen, or stay the same.
///
/// Return if not a leaf; otherwise, apply all Refine refinement
/// criteria to the CommBlock, and set level_desired accordingly:
/// level+1 if it needs to refine, level - 1 if it can coarsen,
/// or level.
///
/// @param[in]  level_maximum   Maximum level to refine
///
/// @return The desired level based on local refinement criteria.
int CommBlock::adapt_compute_desired_level_(int level_maximum)
{

  if (! is_leaf()) return adapt_same;

  adapt_ = adapt_unknown;

  int level = this->level();
  int level_desired = level;

  const FieldDescr * field_descr = simulation()->field_descr();


  Problem * problem = simulation()->problem();
  Refine * refine;

  int index_refine = 0;
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

void CommBlock::adapt_refine_()
{

#ifdef DEBUG_ADAPT
  index_.print("REFINE",-1,2,false,simulation());
#endif

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
  is_leaf_ = false;
#ifdef DEBUG_ADAPT
  index_.print("adapt_refine leaf=0",-1,2,false,simulation());
#endif
}

//----------------------------------------------------------------------

void CommBlock::adapt_delete_child_(Index index_child)
{
#ifdef DEBUG_ADAPT
  index_child.print("adapt_delete_child",-1,2,false,simulation());
#endif
  thisProxy[index_child].p_adapt_delete();

  if (sync_coarsen_.next()) {
    children_.clear();
  }
}

//----------------------------------------------------------------------

void CommBlock::adapt_send_neighbors_levels(int level_new)
{
  if (!is_leaf()) return;

  const int level        = this->level();
  const int rank         = simulation()->dimension();
  const int rank_refresh = simulation()->config()->field_refresh_rank;

  ItFace it_face(rank,rank_refresh);
  int of3[3];

  while (it_face.next(of3)) {

    int ic3[3] = {0,0,0};

    Index index_neighbor = neighbor_(of3);

    // sprintf (buffer,"adapt_send_levels %d %d\n",
    // 	     index_.level(), level);
    // index_.print(buffer,-1,2,false,simulation());
    // index_neighbor.print("neighbor",-1,2,false,simulation());

    const int level_face = face_level (of3);

    if (level_face == level) {

      // SEND-SAME: Face and level are sent to unique
      // neighboring block in the same level

   
      PUT_LEVEL (index_,index_neighbor,ic3,of3,level,level_new,"send-same");

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

	PUT_LEVEL (index_,index_uncle,ic3,of3,level,level_new,"send-coarse");

      }

    } else if (level_face == level + 1) {

      // SEND-FINE: Face and level are sent to all nibling
      // blocks in the next-finer level along the face.

      const int if3[3] = {-of3[0],-of3[1],-of3[2]};
      ItChild it_child(rank,if3);
      while (it_child.next(ic3)) {
	Index index_nibling = index_neighbor.index_child(ic3);

	PUT_LEVEL (index_,index_nibling,ic3,of3,level,level_new,"send-fine");

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

/// @brief Entry function for receiving desired level of a neighbor
///
/// @param index_caller   mesh index of the calling neighbor
/// @param ic3            child indices of neighbor if it's in a finer level
/// @param if3            face (inward) shared with neighbor
/// @param level_face     neighbor's current level
/// @param level_face_new neighbor's desired level
///
/// level_face
/// level_face_new
/// level_new
/// level_new_
/// level
void CommBlock::p_adapt_recv_level
(
 Index index_send,
 int ic3[3],
 int if3[3], 
 int level_face_curr,
 int level_face_new
 )
{

  if (index_send.level() != level_face_curr) {
    PARALLEL_PRINTF ("%d level mismatch between index_send (%d) and level_face_curr (%d)",
		     __LINE__,index_send.level(), level_face_curr);
    index_.print("index_",-1,2,false,simulation());
    index_send.print("index_",-1,2,false,simulation());
  }

#ifdef DEBUG_ADAPT
  std::string bit_str = index_send.bit_string(-1,2);
  sprintf (buffer,"%12s %8s [%d]"
	   "if3 %2d %2d %2d  ic3 %d %d %d  "			
	   "%d <- %d",					
	   "recv",bit_str.c_str(),face_level_last_[IF3(if3)],
	   if3[0],if3[1],if3[2],				
	   ic3[0],ic3[1],ic3[2],				
	   level_face_curr,level_face_new);			
  index_.print(buffer,-1,2,false,simulation());		
#endif

  if (face_level_last_[IF3(if3)] > level_face_new) {
    return;
  } else {
    face_level_last_[IF3(if3)] = level_face_new;
  }

  int level_new = level_new_;

  const int level        = this->level();

  const int of3[3] = {-if3[0],-if3[1],-if3[2] };

  if ( ! is_leaf()) {

    adapt_recv_recurse(if3,ic3,level_face_curr,level_face_new,index_send);

  } else {

    if (level_face_curr == level) {

      adapt_recv_same(of3,level_face_new);

    } else if ( level_face_curr == level + 1 ) {

      adapt_recv_fine(of3,ic3,level_face_new);

    } else if ( level_face_curr == level - 1 ) {

      adapt_recv_coarse(of3,ic3,level_face_new);

    } else  {

      WARNING2 ("CommBlock::notify_neighbor()",
		"level %d and face level %d differ by more than 1",
		level,level_face_curr);
    }

    //
    // If this block wants to coarsen, then 1. all siblings must be able to
    // coarsen as well, and 2. all non-sibling (nephew) must not be going
    // to a level finer than the current level (otherwise a level jump will occur).
    // If either of these cases is true, then change the desired level to the current
    // level (neither coarsen nor refine) and re-send desired level to neighbors

    // The calling block is a sibling if it has the same parent

    const bool is_coarsening = level_new < level;

    const bool is_sibling = (index_.level() > 0 && index_send.level() > 0) ?
      (index_send.index_parent() == index_.index_parent()) : false;

    // The calling block is a nephew if it is a child of a sibling

    const bool is_nephew = (index_.level() > 0 && index_send.level() > 1) ?
      (index_send.index_parent().index_parent() == index_.index_parent()) : false;

    const bool is_finer_neighbor = level_face_new > level_new;

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
#ifdef DEBUG_ADAPT
	sprintf (buffer,"cannot coarsen sibling %d neighbor_finer %d nephew %d",
		 is_sibling,is_finer_neighbor,is_nephew);
	index_.print(buffer,-1,2,false,simulation());
#endif
	level_new = level;
      }
    }

	  
    if (level_new < level_face_new - 1) {

      // restrict new level to within 1 of neighbor
      level_new = level_face_new - 1;

    }

    // notify neighbors if level_new has changed

    if (level_new != level_new_) {
      ASSERT2 ("CommBlock::p_adapt_recv_level()",
	       "level_new %d level_new_ %d\n", level_new,level_new_,
	       level_new > level_new_);
      level_new_ = level_new;
      adapt_send_neighbors_levels(level_new_);
    }

  }
}

//----------------------------------------------------------------------

void CommBlock::adapt_recv_same(const int of3[3],int level_face_new)
{
  // RECV-SAME: Face and level are received from unique
  // neighbor.  Unique face level is updated, and levels on
  // possibly multiple faces of multiple children are updated.

  const int rank = simulation()->dimension();
  const int rank_refresh = simulation()->config()->field_refresh_rank;

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
}

//----------------------------------------------------------------------

void CommBlock::adapt_recv_coarse(const int of3[3], const int ic3[3], int level_face_new)
{      // RECV-COARSE: Face and level are received from unique
  // neighbor.  Possibly multiple faces of block are updated
  // corresponding to the coarse neighbor's face.  Levels of
  // possibly multiple faces of possibly multiple child faces are
  // updated.

  const int rank = simulation()->dimension();
  const int rank_refresh = simulation()->config()->field_refresh_rank;

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
}

//----------------------------------------------------------------------

void CommBlock::adapt_recv_fine(const int of3[3], const int ic3[3],int level_face_new)
{
  // RECV-FINE: Face, level, and sender child indices are received
  // from possibly non-unique neighbor for the corresponding face.
  // If neighbor level indicates the neighbor needs to refine,
  // then refinement is triggered in this block as well.  Face
  // levels are ignored, since the face either stays the same, or
  // refinement is triggered, in which case the child face levels
  // are used.  The neighbor level for the unique face and unique
  // child facing the sending child is updated.

  const int rank = simulation()->dimension();
  const int rank_refresh = simulation()->config()->field_refresh_rank;


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

}

//----------------------------------------------------------------------

void CommBlock::adapt_recv_recurse(const int if3[3], 
				   const int ic3[3], 
				   int level_face_curr, int level_face_new,
				   Index index_send)
{
  const int rank = simulation()->dimension();

#ifdef DEBUG_ADAPT
  sprintf (buffer,"Recurse called: ic3 (%d %d %d) if3 (%d %d %d) level %d -> %d",
	   ic3[0],ic3[1],ic3[2],if3[0],if3[1],if3[2],
	   level_face_curr,level_face_new);
  index_.print(buffer,-1,2,false,simulation());
#endif    

  // Forward to children if internal node
    
  const int of3[3] = {-if3[0],-if3[1],-if3[2] };
  ItChild it_child (rank,of3);

  int jc3[3];
  while (it_child.next(jc3)) {
    Index index_child = index_.index_child(jc3);
    // --------------------------------------------------
    PUT_LEVEL (index_send,index_child,ic3,if3,level_face_curr,level_face_new,"RECURSE");
    // --------------------------------------------------

  }
}

//----------------------------------------------------------------------

void CommBlock::adapt_coarsen_()
{
#ifdef DEBUG_ADAPT
  index_.print("COARSEN",-1,2,false,simulation());
#endif
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
    //    bool lghost[3] = {true,true,true};
    bool lghost[3] = {false,false,false};
    FieldFace * field_face = 
      load_face_(&narray,&array, iface, ic3, lghost, op_array_restrict);

    // copy face levels
    int nf = face_level_.size();
    int face_level[nf];
    for (int i=0; i<nf; i++) face_level[i] = face_level_[i];

    // send child data to parent
    // --------------------------------------------------
    thisProxy[index_parent].p_adapt_recv_child
      (ic3, narray,array, nf,face_level);
    // --------------------------------------------------

    delete field_face;

  }
}

//----------------------------------------------------------------------

void CommBlock::p_adapt_recv_child
(
 int ic3[3],
 int na, char * array,
 int nf, int * child_face_level
 )
{
  // copy array
  int iface[3] = {0,0,0};
  bool lghost[3] = {false,false,false};
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

  is_leaf_=true;
#ifdef DEBUG_ADAPT
  index_.print("p_adapt_recv_child=1",-1,2,false,simulation());
#endif

  adapt_delete_child_(index_child);

  age_ = 0;
}

//----------------------------------------------------------------------

void CommBlock::p_adapt_delete()
{
#ifdef DEBUG_ADAPT
  index_.print("DELETING",-1,2,false,simulation());
#endif
  delete_ = true;
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
  const FieldDescr * field_descr = simulation()->field_descr();
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

