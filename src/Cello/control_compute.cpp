// See LICENSE_CELLO file for license and copyright information

/// @file     control_compute.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ compute-related functions
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//======================================================================

void CommBlock::compute_begin_ ()
{


#ifdef CONFIG_USE_PROJECTIONS
  //  double time_start = CmiWallTimer();
#endif

  if (is_leaf_) {

    const Problem * problem = simulation()->problem();
    Method * method;
    int index_method = 0;
    while ((method = problem->method(index_method++) )) {

      method -> compute (this);

    }
  }

#ifdef CONFIG_USE_PROJECTIONS
  //  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  // Update CommBlock cycle and time to Simulation time and cycle
  // if (cycle_>460) {
  //   char buffer[80];
  //   sprintf (buffer,"time %f dt %f cycle %d\n",time_,dt_,cycle_);
  //   index_.print(buffer,-1,2,false,simulation());
  // }

  set_cycle (cycle_ + 1);
  set_time  (time_  + dt_);
  simulation()->set_cycle(cycle_);
  simulation()->set_time(time_);

  TRACE ("END   PHASE COMPUTE");

  control_sync (sync_compute_exit);
}

//----------------------------------------------------------------------


