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
  CkPrintf ("%d %s:%d %s TRACE_CONTROL %s \n",				\
	    CkMyPe(),__FILE__,__LINE__,					\
	    name_.c_str(), A);						\
  fflush(stdout);						
#else
# define TRACE_CONTROL(A) ;
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

    proxy_simulation[0].p_monitor();
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

  adapt_enter_();
}

//----------------------------------------------------------------------

void Block::refresh_exit_()
{
  TRACE_CONTROL("refresh_exit");

  update_boundary_();

  // CkCallback (refresh_.back()->callback(),thisProxy).send(NULL);

  control_sync (refresh_.back()->callback(),
		refresh_.back()->sync_type(),
		refresh_.back()->sync_exit());
    
#ifdef DEBUG_REFRESH 
 printf ("%d DEBUG_REFRESH Calling Block %s refresh_pop_back(%p)\n",
	  CkMyPe(),name().c_str(),refresh());
  fflush(stdout);
      
#endif
  
  //   Refresh * refresh = refresh_.back();
  //   delete refresh;
  //   refresh_.pop_back();
}

//----------------------------------------------------------------------

void Block::control_sync (int entry_point, int sync, int id)
{
  TRACE_CONTROL("control_sync()");

  if (sync == sync_quiescence) {

    control_sync_quiescence (entry_point);

  } else if (sync == sync_neighbor) {

    control_sync_neighbor (entry_point,id);

  } else if (sync == sync_face) {
 
    control_sync_face (entry_point,id);

  } else if (sync == sync_barrier) {

    control_sync_barrier (entry_point);

  } else {

     ERROR1 ("Block::control_sync()",  "Unknown sync type %d", sync);    

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

void Block::control_sync_neighbor(int entry_point, int id_sync)
{
  TRACE_CONTROL("control_sync_neighbor");

  if ( ! is_leaf() ) {

    TRACE_CONTROL("control_sync_neighbor ! is_leaf()");
    CkCallback(entry_point,CkArrayIndexIndex(index_),thisProxy).send();

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

  const int min_face_rank = 0;

  int num_neighbors = 0;

  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);

  int of3[3];  // ignored
  while (it_neighbor.next(of3)) {

    ++num_neighbors;

    Index index_neighbor = it_neighbor.index();

#ifdef DEBUG_CONTROL
    CkPrintf ("%s calling p_control_sync_count (%d %d 0)\n",
	      name().c_str(),entry_point,id_sync);
    fflush(stdout);
#endif    
    thisProxy[index_neighbor].p_control_sync_count(entry_point,id_sync,0);

  }
#ifdef DEBUG_CONTROL
    CkPrintf ("%s calling p_control_sync_count (%d %d 0)\n",
	      name().c_str(), entry_point,id_sync);
    fflush(stdout);
#endif    
  control_sync_count (entry_point,id_sync,num_neighbors + 1);

}

//----------------------------------------------------------------------

void Block::control_sync_face(int entry_point, int id_sync)
{

  TRACE_CONTROL("control_sync_face");

  const int min_face_rank = 0;

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
  if (id_sync >= sync_max_.size()) {
    sync_count_.resize(id_sync+1);
    sync_max_.resize(id_sync+1);
    sync_count_[id_sync] = 0;
    sync_max_[id_sync] = 0;
  }
#ifdef DEBUG_CONTROL
  CkPrintf ("%s control_sync_count %d %d %d/%d\n",
	    name().c_str(),entry_point,id_sync,count,sync_max_[id_sync]);
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

