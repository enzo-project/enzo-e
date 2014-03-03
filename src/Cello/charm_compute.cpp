// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_compute.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ compute-related functions

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

    FieldDescr * field_descr = simulation()->field_descr();
    const Problem * problem = simulation()->problem();
    Method * method;
    int index_method = 0;
    while ((method = problem->method(index_method++) )) {

      method -> compute_block (field_descr,this);

    }
  }

#ifdef CONFIG_USE_PROJECTIONS
  //  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  // Update CommBlock cycle and time to Simulation time and cycle

  set_cycle (cycle_ + 1);
  set_time  (time_  + dt_);
  
  TRACE ("END   PHASE COMPUTE");

  control_sync (phase_sync_compute_exit);
}

//----------------------------------------------------------------------


