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

  control_sync(CkIndex_Block::r_adapt_enter(NULL),sync_barrier); 
  performance_stop_(perf_initial);
}

//----------------------------------------------------------------------

void Block::adapt_enter_()
{
  TRACE_CONTROL("adapt_enter");
  
  if ( do_adapt_()) {

    adapt_begin_();
    
  } else {

    adapt_exit_();

  }
}

//----------------------------------------------------------------------

void Block::adapt_exit_()
{
  TRACE_CONTROL("adapt_exit");

  control_sync(CkIndex_Main::p_output_enter(),sync_quiescence);
}

//----------------------------------------------------------------------

void Block::output_enter_ ()
{
  performance_start_(perf_output);
  TRACE_CONTROL("output_enter");

  output_begin_();
  performance_stop_(perf_output);
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

  control_sync(CkIndex_Block::r_stopping_enter(NULL),sync_barrier);

}

//----------------------------------------------------------------------

void Block::stopping_enter_()
{

  TRACE_CONTROL("stopping_enter");

  stopping_begin_();
}

//----------------------------------------------------------------------

void Block::stopping_exit_()
{
  TRACE_CONTROL("stopping_exit");

  if (simulation()->cycle_changed()) {
    // if performance counters haven't started yet for this cycle
    int cycle_initial = simulation()->config()->initial_cycle;
    if (cycle_ > cycle_initial) {
      // stop if any previous cycle
      performance_stop_(perf_cycle,__FILE__,__LINE__);
    }
    // start 
    performance_start_ (perf_cycle,__FILE__,__LINE__);
  }

  if (stop_) {

    control_sync(CkIndex_Block::r_exit(NULL),sync_barrier);

  } else {

    compute_enter_();

  }
}

//----------------------------------------------------------------------

void Block::compute_enter_ ()
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  TRACE_CONTROL("compute_enter");

  compute_begin_();
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void Block::compute_exit_ ()
{
  TRACE_CONTROL("compute_exit");

  adapt_enter_();
}

//----------------------------------------------------------------------

void Block::refresh_enter (int callback, Refresh * refresh) 
{
  TRACE_CONTROL("refresh_enter");

#ifdef DEBUG_REFRESH
  CkPrintf ("%s:%d %s Block::set_refresh (%p)\n",
	    __FILE__,__LINE__,name().c_str(),refresh);
  fflush(stdout);
#endif  
  set_refresh(refresh);

  // Update refresh object for the Block

  refresh_.back()->set_callback(callback);

  refresh_begin_();
}

//----------------------------------------------------------------------

void Block::refresh_exit_()
{
  TRACE_CONTROL("refresh_exit");

  update_boundary_();

  control_sync(refresh_.back()->callback(), refresh_.back()->sync_type());
#ifdef DEBUG_REFRESH
  printf ("DEBUG_REFRESH Calling Block %s refresh_pop_back(%p)\n",
	  name().c_str(),refresh());
  fflush(stdout);
      
#endif
  Refresh * refresh = refresh_.back();
  delete refresh;
  refresh_.pop_back();
}

//----------------------------------------------------------------------

void Block::control_sync (int entry_point, int sync, int id)
{
  TRACE_CONTROL("control_sync()");

  if (sync == sync_quiescence) {

    if (index_.is_root())
      CkStartQD(CkCallback (entry_point,proxy_main));

  } else if (sync == sync_neighbor) {

#ifdef DEBUG_REFRESH    
    CkPrintf ("DEBUG_REFRESH %s neighbor sync id %d\n",
	      name().c_str(),id);
#endif    
    control_sync_neighbor_(entry_point,id);

  } else if (sync == sync_face) {
 
    control_sync_face_(entry_point,id);

  } else if (sync == sync_barrier) {

    contribute(CkCallback (entry_point,thisProxy));

  } else {
    ERROR1 ("Block::control_sync()",  "Unknown sync type %d", sync);    
  }
}

//----------------------------------------------------------------------

void Block::control_sync_neighbor_(int entry_point, int phase)
{
  TRACE_CONTROL("control_sync_neighbor");
  
  if ( ! is_leaf() ) {

    TRACE_CONTROL("control_sync_neighbor ! is_leaf()");
    CkCallback(entry_point,CkArrayIndexIndex(index_),thisProxy).send();

    return;
  }

  const int min_face_rank = 0;

  int num_neighbors = 0;

  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);

  while (it_neighbor.next()) {

    ++num_neighbors;

    Index index_neighbor = it_neighbor.index();

#ifdef DEBUG_CONTROL
    CkPrintf ("%s calling p_control_sync_count (%d %d 0)\n",
	      name().c_str(),entry_point,phase);
    fflush(stdout);
#endif    
    thisProxy[index_neighbor].p_control_sync_count(entry_point,phase,0);

  }
#ifdef DEBUG_CONTROL
    CkPrintf ("%s calling p_control_sync_count (%d %d 0)\n",
	      name().c_str(), entry_point,phase);
    fflush(stdout);
#endif    
  control_sync_count_(entry_point,phase,num_neighbors + 1);

}

//----------------------------------------------------------------------

void Block::control_sync_face_(int entry_point, int phase)
{

  TRACE_CONTROL("control_sync_face");

  const int min_face_rank = 0;

  int num_faces = 0;

  ItFace it_face = this->it_face(min_face_rank,index_);

  while (it_face.next()) {

    int of3[3];
    it_face.face(of3);
    
    // Only count face if a Block exists in the level
    if (face_level(of3) >= level()) {
      ++num_faces;

      Index index_face = it_face.index();

      thisProxy[index_face].p_control_sync_count(entry_point,phase,0);
    }

  }
  control_sync_count_(entry_point,phase,num_faces + 1);
}

//----------------------------------------------------------------------

void Block::control_sync_count_ (int entry_point, int phase, int count)
{
  if (phase >= sync_max_.size()) {
    sync_count_.resize(phase+1);
    sync_max_.resize(phase+1);
    sync_count_[phase] = 0;
    sync_max_[phase] = 0;
  }
#ifdef DEBUG_CONTROL
  CkPrintf ("%s control_sync_count %d %d %d/%d\n",
	    name().c_str(),entry_point,phase,count,sync_max_[phase]);
  fflush(stdout);
#endif
  
  if (count != 0)  sync_max_[phase] = count;

  ++sync_count_[phase];

  // sync_max_ reached: continue and reset counter

  if (sync_max_[phase] > 0 && sync_count_[phase] >= sync_max_[phase]) {

    sync_max_  [phase] = 0;
    sync_count_[phase] = 0;

    CkCallback(entry_point,CkArrayIndexIndex(index_),thisProxy).send(NULL);

  }
}

//======================================================================

