// See LICENSE_CELLO file for license and copyright information

/// @file     control_charm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-02-13
/// @brief    Functions controling control flow of charm entry functions
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_REFRESH
// #define DEBUG_CONTROL

#ifdef DEBUG_CONTROL
# define TRACE_CONTROL(A)						\
  CkPrintf ("%d %s %s TRACE_CONTROL %s \n",				\
	    CkMyPe(),__FILE__,					\
	    name_.c_str(), A);						\
  fflush(stdout);						
# define TRACE_SYNC(A)							\
  CkPrintf ("%d %s %s TRACE_CONTROL %s entry %d id %d\n",	\
	    CkMyPe(),__FILE__,					\
	    name_.c_str(), A,entry_point,id_sync);			\
  fflush(stdout);						
#else
# define TRACE_CONTROL(A) ;
# define TRACE_SYNC(A) ;
#endif

//----------------------------------------------------------------------

void Block::initial_exit_()
{
  performance_start_(perf_initial);
  TRACE_CONTROL("initial_exit");

#ifdef TRACE_CONTRIBUTE  
  CkPrintf ("%s %s:%d DEBUG_CONTRIBUTE calling r_adapt_enter\n",
	    name().c_str(),__FILE__,__LINE__); fflush(stdout);
#endif  
  control_sync_barrier (CkIndex_Block::r_adapt_enter(NULL));
  performance_stop_(perf_initial);
}

//----------------------------------------------------------------------

void Block::adapt_exit_()
{
  TRACE_CONTROL("adapt_exit");

  control_sync_quiescence(CkIndex_Main::p_output_enter());
}

//----------------------------------------------------------------------

void Block::output_exit_()
{
  performance_start_(perf_output);

  TRACE_CONTROL("output_exit");

  if (index_.is_root()) {
    cello::simulation()->monitor_output();
  }

  performance_stop_(perf_output);

#ifdef TRACE_CONTRIBUTE  
  CkPrintf ("%s %s:%d DEBUG_CONTRIBUTE calling r_stopping_enter()\n",
	    name().c_str(),__FILE__,__LINE__); fflush(stdout);
#endif  
  control_sync_barrier (CkIndex_Block::r_stopping_enter(NULL));

}

//----------------------------------------------------------------------

void Block::stopping_exit_()
{
  TRACE_CONTROL("stopping_exit");

  if (cello::simulation()->cycle_changed()) {
    // if performance counters haven't started yet for this cycle
    int cycle_initial = cello::config()->initial_cycle;
    if (cycle_ > cycle_initial) {
      // stop if any previous cycle
      performance_stop_(perf_cycle,__FILE__,__LINE__);
    }
    // start 
    performance_start_ (perf_cycle,__FILE__,__LINE__);
  }

  if (stop_) {

#ifdef TRACE_CONTRIBUTE  
  CkPrintf ("%s %s:%d DEBUG_CONTRIBUTE calling r_exit()\n",
	    name().c_str(),__FILE__,__LINE__); fflush(stdout);
#endif  
    control_sync_barrier (CkIndex_Block::r_exit(NULL));

  } else {

    compute_enter_();

  }
}

//----------------------------------------------------------------------

void Block::compute_exit_ ()
{
  TRACE_CONTROL("compute_exit");

  control_sync_barrier(CkIndex_Block::r_adapt_enter(NULL));
  //  adapt_enter_();
}

//----------------------------------------------------------------------

void Block::refresh_exit_()
{
  TRACE_CONTROL("refresh_exit");

  update_boundary_();

  Refresh * refresh = refresh_.back();

#ifdef NEW_REFRESH
  // control_sync (refresh->callback(),
  // 		refresh->sync_type(),
  // 		refresh->sync_exit(),
  // 		refresh->min_face_rank(),
  // 		refresh->neighbor_type(),
  // 		refresh->root_level());
  ERROR1 ("Block::refresh_exit_",
	  "%s Should not be called with NEW_REFRESH",this->name().c_str());
  
  if (index().is_root()) {
    WARNING("Block::refresh_exit_()",
  	    "remove sync from refresh_exit_() callback with NEW_REFRESH");
  }
  CkCallback
    (refresh->callback(),
     CkArrayIndexIndex(index_),thisProxy).send(NULL);
#else
#ifdef DEBUG_REFRESH 
  printf ("%d DEBUG_REFRESH refresh_exit calling refresh\n",
	  CkMyPe(),name().c_str());
  fflush(stdout);
      
#endif

  control_sync (refresh->callback(),
  		refresh->sync_type(),
  		refresh->sync_exit(),
  		refresh->min_face_rank(),
  		refresh->neighbor_type(),
  		refresh->root_level());
#endif  
    
  // WARNING: BREAKS Charm++ with random queueing on some regression tests
  //  delete refresh;
  //  refresh_.pop_back();
}

//----------------------------------------------------------------------

void Block::control_sync (int entry_point, int sync_type, int id_sync,
			  int min_face_rank, int neighbor_type, int root_level)
{
  TRACE_CONTROL("control_sync()");
  TRACE_SYNC("control_sync()");
#ifdef DEBUG_CONTROL
  CkPrintf ("DEBUG_CONTROL %s entry_point = %d\n",name().c_str(),entry_point);
  fflush(stdout);
#endif  

  if (sync_type == sync_quiescence) {

    control_sync_quiescence (entry_point);

  } else if (sync_type == sync_neighbor) {

    control_sync_neighbor
      (entry_point,id_sync,min_face_rank, neighbor_type,root_level);

  } else if (sync_type == sync_face) {
 
    control_sync_face (entry_point,id_sync,min_face_rank);

  } else if (sync_type == sync_barrier) {

    control_sync_barrier (entry_point);

  } else {

     ERROR2 ("Block::control_sync()",
	     "Unknown sync type %d neighbor type %d",
	     sync_type,neighbor_type);    

  }
}

//----------------------------------------------------------------------

void Block::control_sync_quiescence (int entry_point)
{
  if (index_.is_root())
    CkStartQD(CkCallback (entry_point,proxy_main));
}

//----------------------------------------------------------------------

void Block::control_sync_barrier (int entry_point)
{
  contribute(CkCallback (entry_point,thisProxy));
}

//----------------------------------------------------------------------

void Block::control_sync_neighbor(int entry_point, int id_sync,
				  int min_face_rank,
				  int neighbor_type,
				  int root_level)
{
  TRACE_CONTROL("control_sync_neighbor");
  TRACE_SYNC("control_sync_neighbhor()");

  if ( ! is_leaf() ) {

    TRACE_CONTROL("control_sync_neighbor ! is_leaf()");
    CkCallback(entry_point,CkArrayIndexIndex(index_),thisProxy).send(NULL);

    return;
  }

  ASSERT1("control_sync()",
	  "id %d must be specified for neighbor sync and >= 0",
	  id_sync,
	  (id_sync >= 0));

#ifdef DEBUG_REFRESH    
  CkPrintf ("%d DEBUG_REFRESH %s neighbor sync id %d\n",
	    CkMyPe(), name().c_str(),id_sync);
#endif    

  int num_neighbors = 0;

  const int min_level = cello::config()->mesh_min_level;
  
  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_,neighbor_type,min_level,root_level);

  int of3[3];  // ignored
  while (it_neighbor.next(of3)) {

    ++num_neighbors;

    Index index_neighbor = it_neighbor.index();

#ifdef DEBUG_CONTROL
    CkPrintf ("%s DEBUG_CONTROL calling p_control_sync_count (%d %d 0)\n",
	      name().c_str(),entry_point,id_sync);
    fflush(stdout);
#endif    
    thisProxy[index_neighbor].p_control_sync_count(entry_point,id_sync,0);

  }
#ifdef DEBUG_CONTROL
    CkPrintf ("%s DEBUG_CONTROL calling p_control_sync_count count %d (%d %d 0)\n",
	      name().c_str(),num_neighbors, entry_point,id_sync);
    fflush(stdout);
#endif    
  control_sync_count (entry_point,id_sync,num_neighbors + 1);

}

//----------------------------------------------------------------------

void Block::control_sync_face(int entry_point, int id_sync, int min_face_rank)
{

  TRACE_CONTROL("control_sync_face");
  TRACE_SYNC("control_sync_face()");

  int num_faces = 0;

  ItFace it_face = this->it_face(min_face_rank,index_);

  int of3[3];
  while (it_face.next(of3)) {

    // Only count face if a Block exists in the level
    if (face_level(of3) >= level()) {
      ++num_faces;

      Index index_face = it_face.index();

      thisProxy[index_face].p_control_sync_count(entry_point,id_sync,0);
    }

  }
  control_sync_count (entry_point,id_sync,num_faces + 1);
}

//----------------------------------------------------------------------

void Block::control_sync_count (int entry_point, int id_sync, int count)
{
  const int n_new = id_sync + 1;
  const int n_now = sync_max_.size();
  if (n_new  > n_now) {
    sync_count_.resize(n_new,0);
    sync_max_.resize(n_new,0);
  }
#ifdef DEBUG_CONTROL
  CkPrintf ("%s DEBUG_CONTROL control_sync_count %d %d %d/%d\n",
	    name().c_str(),entry_point,id_sync,sync_count_[id_sync],sync_max_[id_sync]);
  fflush(stdout);
#endif
  
  if (count != 0)  sync_max_[id_sync] = count;

  ++sync_count_[id_sync];

  // sync_max_ reached: continue and reset counter

  if (sync_max_[id_sync] > 0 && sync_count_[id_sync] >= sync_max_[id_sync]) {

    sync_max_  [id_sync] = 0;
    sync_count_[id_sync] = 0;

    CkCallback(entry_point,CkArrayIndexIndex(index_),thisProxy).send(NULL);

  }
}

//======================================================================

