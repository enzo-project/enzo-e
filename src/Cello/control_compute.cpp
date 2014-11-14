// See LICENSE_CELLO file for license and copyright information

/// @file     control_compute.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ compute-related functions
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//======================================================================

void CommBlock::compute_begin_ ()
{
  // REFRESH HERE WITH control_

  refresh_phase_ = phase_compute_continue;
  refresh_sync_  = "contribute";

  control_next(phase_refresh_enter,"neighbor");

}

//----------------------------------------------------------------------

void CommBlock::compute_continue_ ()
{

#ifdef CONFIG_USE_PROJECTIONS
  //  double time_start = CmiWallTimer();
#endif

  simulation()->set_phase(phase_compute);

  const Problem * problem = simulation()->problem();
  int index_method = 0;
  Method * method;
  while ((method = problem->method(index_method++) )) {

    // Remember currently active Method to returne from reductions
    set_method (method);

    // Apply the method to the CommBlock
    method -> compute (this);

  }

#ifdef CONFIG_USE_PROJECTIONS
  //  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  set_cycle (cycle_ + 1);
  set_time  (time_  + dt_);
  simulation()->set_cycle(cycle_);
  simulation()->set_time(time_);

  TRACE ("END   PHASE COMPUTE");

  //  control_sync (phase_compute_exit,"none",true,__FILE__,__LINE__);
  control_next();

}

//----------------------------------------------------------------------


