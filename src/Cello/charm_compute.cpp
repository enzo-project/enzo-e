// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_compute.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ compute-related functions
///
/// COMPUTE
///
///    SimulationCharm::c_compute()
///       if (!stop)
///          if (cycle > 0) performance.stop_region()
///          performance.start_region()
///          if (is_root())
///             block_array().p_compute()
///       else
///          performance.stop_region()
///          performance_write()
///	  Main::p_exit()
///
///    CommBlock::p_compute()
///       compute()
///
///    CommBlock::compute()
///       for method
///           method->compute_block()
///       set_state(cycle+1, time + dt)
///
///       if (adapt)
///          >>>>> adapt(0) >>>>>
///       else
///          >>>>> refresh() >>>>> 

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//======================================================================


void SimulationCharm::c_compute()
{
  TRACE("SimulationCharm::c_compute()");
  if (stop_) {
    
    performance_.stop_region (id_cycle_);
    performance_write();

    proxy_main.p_exit(CkNumPes());

  } else {

    if (cycle_ > 0 ) performance_.stop_region (id_cycle_);

    performance_.start_region (id_cycle_);

    if (hierarchy()->group_process()->is_root()) 
      hierarchy()->block_array()->p_compute(cycle_,time_,dt_);
  }
}

//----------------------------------------------------------------------

void CommBlock::p_compute (int cycle, double time, double dt)
{
  // set_cycle(cycle);
  // set_time(time);
  // set_dt(dt);

  TRACE3 ("CommBlock::p_compute() cycle %d time %f dt %f",cycle,time,dt);
  compute();
}

//----------------------------------------------------------------------

void CommBlock::compute()
{
  TRACE ("CommBlock::compute()");

  Simulation * simulation = proxy_simulation.ckLocalBranch();

 #ifdef CONFIG_USE_PROJECTIONS
   double time_start = CmiWallTimer();
 #endif

  FieldDescr * field_descr = simulation->field_descr();

  int index_method = 0;
  while (Method * method = simulation->problem()->method(index_method++)) {
    method -> compute_block (field_descr,this);
  }

 #ifdef CONFIG_USE_PROJECTIONS
   traceUserBracketEvent(10,time_start, CmiWallTimer());
 #endif

  // Update CommBlock cycle and time to Simulation time and cycle

  set_cycle (cycle_ + 1);

  set_time  (time_  + dt_);
  
  // prepare for next cycle: Timestep, Stopping, Monitor, Output

  TRACE ("CommBlock::compute() calling p_adapt()");
  p_adapt(0);

}

//----------------------------------------------------------------------

#endif /* CONFIG_USE_CHARM */


