// See LICENSE_CELLO file for license and copyright information

/// @file     charm_stopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
///
///    STOPPING
///        
///    CommBlock::stopping()
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

//----------------------------------------------------------------------

void CommBlock::stopping_enter_()
{
#ifdef TRACE_MEMORY
  trace_mem_ = Memory::instance()->bytes();
#endif

  performance_switch_(perf_stopping,__FILE__,__LINE__);

  TRACE1("CommBlock::stopping() %p",&thisProxy);
  Simulation * simulation = proxy_simulation.ckLocalBranch();

 //--------------------------------------------------
  // Enforce boundary conditions
  //--------------------------------------------------

  update_boundary_();

  int stopping_interval = simulation->config()->stopping_interval;
  bool stopping_reduce = stopping_interval ? 
    ((cycle_ % stopping_interval) == 0) : false;

  if (stopping_reduce || dt_==0.0) {
    // Compute local dt

    Problem * problem = simulation->problem();

    Timestep * timestep = problem->timestep();

    const FieldDescr * field_descr = simulation->field_descr();
    double dt_block = timestep->evaluate(field_descr,this);
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

    // Evaluate local stopping criteria

    int stop_block = stopping->complete(cycle_,time_);

    // Reduce to find CommBlock array minimum dt and stopping criteria

    double min_reduce[2];

    min_reduce[0] = dt_block;
    min_reduce[1] = stop_block ? 1.0 : 0.0;

    // --------------------------------------------------
    // ENTRY: #1 CommBlock::stopping()-> CommBlock::r_stopping_compute_timestep()
    // ENTRY: contribute()
    // --------------------------------------------------
    CkCallback callback (CkIndex_CommBlock::r_stopping_compute_timestep(NULL), thisProxy);
    contribute(2*sizeof(double), min_reduce, CkReduction::min_double, callback);
    // --------------------------------------------------

  } else {

    stopping_exit_();

  }

}

//----------------------------------------------------------------------

void CommBlock::r_stopping_compute_timestep(CkReductionMsg * msg)
{
  double * min_reduce = (double * )msg->getData();

  dt_   = min_reduce[0];
  stop_ = min_reduce[1] == 1.0 ? true : false;

  delete msg;

  stopping_exit_();
}

//----------------------------------------------------------------------

void CommBlock::stopping_exit_()
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  PARALLEL_PRINTF ("%d dt = %d\n",CkMyPe(),dt_);
  set_dt   (dt_);
  set_stop (stop_);

  simulation->update_state(cycle_,time_,dt_,stop_);

  SimulationCharm * simulation_charm = proxy_simulation.ckLocalBranch();

#ifdef TRACE_MEMORY
  trace_mem_ = Memory::instance()->bytes() - trace_mem_;
  PARALLEL_PRINTF ("memory stopping %lld\n",trace_mem_);
#endif

  simulation_charm->begin_output();
}

