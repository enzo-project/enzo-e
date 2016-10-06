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

//======================================================================

void Block::compute_begin_ ()
{
#ifdef DEBUG_COMPUTE
  CkPrintf ("%s DEBUG_COMPUTE Block::compute_begin_()\n", name().c_str());
#endif

  simulation()->set_phase(phase_compute);

  index_method_ = 0;
  compute_next_();
}

//----------------------------------------------------------------------

void Block::compute_next_ ()
{
#ifdef DEBUG_COMPUTE
  CkPrintf ("%s DEBUG_COMPUTE Block::compute_next_()\n", name().c_str());
#endif

  Method * method = this->method();

  if (method) {

    Refresh * refresh = method->refresh();

    refresh->set_active (is_leaf());

    refresh_enter
      (CkIndex_Block::r_compute_continue(NULL), refresh );

  } else {

    compute_end_();

  }
}

//----------------------------------------------------------------------

void Block::compute_continue_ ()
{
#ifdef DEBUG_COMPUTE
  CkPrintf ("%s DEBUG_COMPUTE Block::compute_continue_()\n", name().c_str());
#endif

#ifdef CONFIG_USE_PROJECTIONS
  //  double time_start = CmiWallTimer();
#endif

    Method * method = this->method();
    Schedule * schedule = method->schedule();
    bool is_scheduled = 
      (!schedule) ||
      (schedule && schedule->write_this_cycle(cycle_,time_));

    // printf ("DEBUG Method %s schedule = %p scheduled %d\n",
    // 	    method->name().c_str(),
    // 	    schedule,
    // 	    is_scheduled);

    if (is_scheduled) {


    TRACE2 ("Block::compute_continue() method = %d %p\n",
	    index_method_,method); fflush(stdout);

#ifdef DEBUG_COMPUTE
    CkPrintf ("%s DEBUG_COMPUTE applying Method %s\n",
	      name().c_str(),method->name().c_str());
#endif
    // Apply the method to the Block

    method -> compute (this);

  } else {

    compute_done();

  }
}

//----------------------------------------------------------------------

void Block::compute_done ()
{
#ifdef DEBUG_COMPUTE
  CkPrintf ("%s DEBUG_COMPUTE Block::compute_done_()\n", name().c_str());
#endif
  index_method_++;
  compute_next_();
}

//----------------------------------------------------------------------

void Block::compute_end_ ()
{
#ifdef DEBUG_COMPUTE
  CkPrintf ("%s DEBUG_COMPUTE Block::compute_end_()\n", name().c_str());
#endif


#ifdef CONFIG_USE_PROJECTIONS
  //  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  set_cycle (cycle_ + 1);
  set_time  (time_  + dt_);
  simulation()->set_cycle(cycle_);
  simulation()->set_time(time_);

  compute_exit_();

  TRACE ("END   PHASE COMPUTE");
}

//----------------------------------------------------------------------



