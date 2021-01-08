// See LICENSE_CELLO file for license and copyright information

/// @file     control_compute.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ compute-related functions
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_COMPUTE

#define CYCLE 0

//======================================================================

void Block::compute_enter_ ()
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  compute_begin_();
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void Block::compute_begin_ ()
{
#ifdef DEBUG_COMPUTE
  if (cycle() >= CYCLE)
    CkPrintf ("%d %s DEBUG_COMPUTE Block::compute_begin_()\n", CkMyPe(),name().c_str());
#endif

  cello::simulation()->set_phase(phase_compute);

  index_method_ = 0;
  compute_next_();
}

//----------------------------------------------------------------------

void Block::compute_next_ ()
{
  Method * method = this->method();

#ifdef DEBUG_COMPUTE
  if (cycle() >= CYCLE)
    CkPrintf ("%d %s DEBUG_COMPUTE Block::compute_next_(%s)\n",CkMyPe(), name().c_str(),method?method->name().c_str():"NULL");
#endif

  if (method) {

#ifdef DEBUG_COMPUTE
    CkPrintf ("DEBUG_REFRESH %s:%d calling refresh_[enter|start]\n",__FILE__,__LINE__);
#endif

    int ir_post = method->refresh_id_post();
    
    cello::refresh(ir_post)->set_active (is_leaf());
    
    new_refresh_start (ir_post,CkIndex_Block::p_compute_continue());
    
  } else {

    compute_end_();

  }
}

//----------------------------------------------------------------------

void Block::compute_continue_ ()
{
  //  this->debug_new_refresh(__FILE__,__LINE__);
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_COMPUTE
  if (cycle() >= CYCLE)
    CkPrintf ("%d %s DEBUG_COMPUTE Block::compute_continue_()\n", CkMyPe(),name().c_str());
#endif

#ifdef CONFIG_USE_PROJECTIONS
  //  double time_start = CmiWallTimer();
#endif

  Method * method = this->method();
  Schedule * schedule = method->schedule();
  bool is_scheduled = 
    (schedule==NULL) ||
    (schedule->write_this_cycle(cycle_,time_));

  if (is_scheduled) {

    TRACE2 ("Block::compute_continue() method = %d %p\n",
	    index_method_,method); fflush(stdout);

#ifdef DEBUG_COMPUTE
    if (cycle() >= CYCLE)
      CkPrintf ("%d %s DEBUG_COMPUTE applying Method %s\n",
	      CkMyPe(),name().c_str(),method->name().c_str());
    CkPrintf ("DEBUG_TRACE_REFRESH Method %s compute()\n",method->name().c_str());
#endif
    // Apply the method to the Block

    method->compute (this);
    
    performance_stop_(perf_compute,__FILE__,__LINE__);

  } else {

    performance_stop_(perf_compute,__FILE__,__LINE__);
    compute_done();

  }
}

//----------------------------------------------------------------------

void Block::compute_done ()
{
#ifdef DEBUG_COMPUTE
  if (cycle() >= CYCLE)
    CkPrintf ("%d %s DEBUG_COMPUTE Block::compute_done_()\n", CkMyPe(),name().c_str());
#endif
  index_method_++;
  compute_next_();
}

//----------------------------------------------------------------------

void Block::compute_end_ ()
{
#ifdef DEBUG_COMPUTE
  if (cycle() >= CYCLE)
    CkPrintf ("%d %s DEBUG_COMPUTE Block::compute_end_()\n", CkMyPe(),name().c_str());
#endif


#ifdef CONFIG_USE_PROJECTIONS
  //  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  // Push back fields if saving old ones
  data()->field().save_history(time_);

  // Update block cycle and time
  set_cycle (cycle_ + 1);
  set_time  (time_  + dt_);

  // Update Simulation cycle and time (redundant)
  cello::simulation()->set_cycle(cycle_);
  cello::simulation()->set_time(time_);

  compute_exit_();

  TRACE ("END   PHASE COMPUTE");
}

//----------------------------------------------------------------------



