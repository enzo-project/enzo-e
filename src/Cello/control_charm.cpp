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

// #define DEBUG_CONTROL

#ifdef DEBUG_CONTROL
# define TRACE_CONTROL(A)						\
  CkPrintf ("%d %s:%d %s TRACE %s\n",					\
	    CkMyPe(),__FILE__,__LINE__,name_.c_str(),A);		\
  fflush(stdout);						
#else
# define TRACE_CONTROL(A) ;
#endif

//----------------------------------------------------------------------

void Block::initial_exit_()
{
  TRACE_CONTROL("initial_exit");

  control_sync(CkIndex_Block::r_adapt_enter(NULL),sync_barrier); 
}

//----------------------------------------------------------------------

void Block::adapt_enter_()
{
  TRACE_CONTROL("adapt_enter");

  performance_switch_ (perf_adapt,__FILE__,__LINE__);

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
  TRACE_CONTROL("output_enter");

  performance_switch_ (perf_output,__FILE__,__LINE__);

  output_begin_();
}

//----------------------------------------------------------------------

void Block::output_exit_()
{

  TRACE_CONTROL("output_exit");

  
  if (index_.is_root()) {

    proxy_simulation[0].p_monitor();
  }

  control_sync(CkIndex_Block::r_stopping_enter(NULL),sync_barrier);
}

//----------------------------------------------------------------------

void Block::stopping_enter_()
{

  TRACE_CONTROL("stopping_enter");

  performance_switch_(perf_stopping,__FILE__,__LINE__);

  stopping_begin_();

}

//----------------------------------------------------------------------

void Block::stopping_exit_()
{

  TRACE_CONTROL("stopping_exit");

  if (stop_) {

    control_sync(CkIndex_Block::r_exit(NULL),sync_barrier);

  } else {

    compute_enter_();

  }

}

//----------------------------------------------------------------------

void Block::compute_enter_ ()
{
  TRACE_CONTROL("compute_enter");

  performance_switch_(perf_compute,__FILE__,__LINE__);

  compute_begin_();
}

//----------------------------------------------------------------------

void Block::compute_exit_ ()
{
  TRACE_CONTROL("compute_exit");

  adapt_enter_();
}

//----------------------------------------------------------------------

void Block::refresh_enter_(int callback, Refresh * refresh) 
{
  TRACE_CONTROL("refresh_enter");
  performance_switch_(perf_refresh,__FILE__,__LINE__);

  set_refresh(refresh);

  // Update refresh object for the Block

  refresh_.set_callback(callback);

  refresh_begin_();
}

//----------------------------------------------------------------------

void Block::refresh_exit_()
{
  TRACE_CONTROL("refresh_exit");

  update_boundary_();

  control_sync(refresh_.callback(), refresh_.sync_type());
}

//----------------------------------------------------------------------

void Block::control_sync (int entry_point, int sync, int id)
{
  
  if (sync == sync_quiescence) {

    if (index_.is_root())
      CkStartQD(CkCallback (entry_point,proxy_main));

  } else if (sync == sync_neighbor) {
 
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
  if ( ! is_leaf() ) {

    CkCallback(entry_point,CkArrayIndexIndex(index_),thisProxy).send();

    return;
  }

  const int min_face_rank = 0;

  int num_neighbors = 0;

  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);

  while (it_neighbor.next()) {

    ++num_neighbors;

    Index index_neighbor = it_neighbor.index();

    thisProxy[index_neighbor].p_control_sync_count(entry_point,phase);

  }
  control_sync_count_(entry_point,phase,num_neighbors + 1);

}

//----------------------------------------------------------------------

void Block::control_sync_face_(int entry_point, int phase)
{

  const int min_face_rank = 0;

  ItFace it_face = this->it_face(min_face_rank,index_);

  int num_faces = 0;
  while (it_face.next()) {

    ++num_faces;

    Index index_face = it_face.index();

    thisProxy[index_face].p_control_sync_count(entry_point,phase);

  }
  control_sync_count_(entry_point,phase,num_faces + 1);

}

//----------------------------------------------------------------------

void Block::control_sync_count_ (int entry_point, int phase, int count)
{
  TRACE_CONTROL("control_sync_count");

  if (count != 0)  max_sync_[phase] = count;

  ++count_sync_[phase];

  // max_sync reached: continue and reset counter

  if (max_sync_[phase] > 0 && count_sync_[phase] >= max_sync_[phase]) {

    max_sync_[phase] = 0;

    count_sync_[phase] = 0;

    CkCallback(entry_point,CkArrayIndexIndex(index_),thisProxy).send();

  }
}

//======================================================================

