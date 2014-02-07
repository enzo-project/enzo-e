// See LICENSE_CELLO file for license and copyright information

/// @file     charm_prepare.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
///
///    PREPARE
///        
///    CommBlock::prepare()
///       update_boundary_()
///       compute dt
///       compute stopping
///       contribute( >>>>> CommBlock::r_output() >>>>> )

#define TRACE_CELLO
   
#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

void CommBlock::prepare()
{
  // index_.print("BEGIN prepare()",-1,2,false,simulation());
  switch_performance_(perf_prepare,__FILE__,__LINE__);

  TRACE1("CommBlock::prepare() %p",&thisProxy);
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  
 //--------------------------------------------------
  // Enforce boundary conditions
  //--------------------------------------------------

  update_boundary_();

  FieldDescr * field_descr = simulation->field_descr();

  //--------------------------------------------------
  // Compute local dt
  //--------------------------------------------------

  Problem * problem = simulation->problem();

  double dt_block;
  Timestep * timestep = problem->timestep();

  dt_block = timestep->evaluate(field_descr,this);
  TRACE1("dt_block = %f",dt_block);

  // Reduce timestep to coincide with scheduled output if needed

  int index_output=0;
  while (Output * output = problem->output(index_output++)) {
    Schedule * schedule = output->schedule();
    dt_block = schedule->update_timestep(time_,dt_block);
  }


  // Reduce timestep to not overshoot final time from stopping criteria

  Stopping * stopping = problem->stopping();

  double time_stop = stopping->stop_time();
  double time_curr = time_;

  dt_block = MIN (dt_block, (time_stop - time_curr));

  //--------------------------------------------------
  // Evaluate local stopping criteria
  //--------------------------------------------------

  int stop_block = stopping->complete(cycle_,time_);

  //--------------------------------------------------
  // Reduce to find CommBlock array minimum dt and stopping criteria
  //--------------------------------------------------

  double min_reduce[2];

  min_reduce[0] = dt_block;
  min_reduce[1] = stop_block ? 1.0 : 0.0;
  
  CkCallback callback (CkIndex_CommBlock::r_output(NULL), thisProxy);

  // --------------------------------------------------
  // ENTRY: #1 CommBlock::prepare()-> CommBlock::r_output()
  // ENTRY: contribute()
  // --------------------------------------------------
  contribute( 2*sizeof(double), min_reduce, CkReduction::min_double, callback);
  // --------------------------------------------------

  //  performance->stop_region(perf_prepare);
}
