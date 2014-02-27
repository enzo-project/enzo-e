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
//----------------------------------------------------------------------

void CommBlock::compute_enter_ (int cycle, double time, double dt)
{
// #ifdef CELLO_TRACE
//   index_.print("BEGIN PHASE COMPUTE p_compute_enter()",-1,2,false,simulation());
// #endif

  // set_cycle(cycle);
  // set_time(time);
  // set_dt(dt);

  performance_switch_(perf_compute,__FILE__,__LINE__);

  TRACE3 ("CommBlock::p_compute_enter() cycle %d time %f dt %f",cycle,time,dt);

#ifdef CONFIG_USE_PROJECTIONS
  //  double time_start = CmiWallTimer();
#endif

  if (is_leaf_) {

// #ifdef CELLO_TRACE
//     index_.print("p_compute_enter",-1,2,false,simulation());
// #endif    
    FieldDescr * field_descr = simulation()->field_descr();
    int index_method = 0;
    while (Method * method = simulation()->problem()->method(index_method++)) {

      method -> compute_block (field_descr,this);

    }
  }

#ifdef CONFIG_USE_PROJECTIONS
  //  traceUserBracketEvent(10,time_start, CmiWallTimer());
#endif

  // Update CommBlock cycle and time to Simulation time and cycle

  set_cycle (cycle_ + 1);
  TRACE1("Calling set_time (%f)",time_+dt_);
  set_time  (time_  + dt_);
  
  TRACE ("CommBlock::compute() calling p_adapt(0)");

  TRACE ("END   PHASE COMPUTE");

  int adapt_interval = simulation()->config()->mesh_adapt_interval;
  if (adapt_interval && ((cycle_ % adapt_interval) == 0)) {
    next_phase_ = phase_adapt;
  } else {
    next_phase_ = phase_stopping;
  }

  compute_exit_();
}

//----------------------------------------------------------------------

void CommBlock::compute_exit_ ()
{
  refresh_enter_();
}

//----------------------------------------------------------------------


