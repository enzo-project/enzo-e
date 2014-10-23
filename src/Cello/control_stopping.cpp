// See LICENSE_CELLO file for license and copyright information

/// @file     control_stopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
/// @ingroup  Control
///
///    STOPPING
///        
///    CommBlock::stopping()
///       update_boundary_()
///       compute dt
///       compute stopping
///       contribute( >>>>> CommBlock::r_output() >>>>> )

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define CELLO_VERBOSE

#ifdef CELLO_VERBOSE
#   define VERBOSE(A)					\
  if (index_.is_root()) {				\
    Monitor * monitor = simulation()->monitor();	\
    monitor->print("Control", A);			\
  } 
#else
#   define VERBOSE(A) ;
#endif


//----------------------------------------------------------------------

void CommBlock::stopping_begin_()
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  simulation->set_phase(phase_stopping);

  int stopping_interval = simulation->config()->stopping_interval;

  bool stopping_reduce = stopping_interval ? 
    ((cycle_ % stopping_interval) == 0) : false;

  if (stopping_reduce || dt_==0.0) {

    // Compute local dt

    Problem * problem = simulation->problem();

    int index = 0;
    Method * method;
    double dt_block = std::numeric_limits<double>::max();
    while ((method = problem->method(index++))) {
      dt_block = std::min(dt_block,method->timestep(this));
    }

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

    CkCallback callback (CkIndex_CommBlock::r_stopping_compute_timestep(NULL),
			 thisProxy);

    contribute(2*sizeof(double), min_reduce, CkReduction::min_double, callback);

  } else {

    stopping_balance_();

  }

}

//----------------------------------------------------------------------

void CommBlock::r_stopping_compute_timestep(CkReductionMsg * msg)
{

  ++age_;

  double * min_reduce = (double * )msg->getData();

  dt_   = min_reduce[0];
  stop_ = min_reduce[1] == 1.0 ? true : false;

  delete msg;

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  set_dt   (dt_);
  set_stop (stop_);

  simulation->set_dt(dt_);
  simulation->set_stop(stop_);

  if (cycle_ > 0 ) {
    performance_stop_(perf_cycle,__FILE__,__LINE__);
  }
  performance_start_ (perf_cycle,__FILE__,__LINE__);

#ifdef CONFIG_USE_PROJECTIONS  
  bool was_off = (simulation->projections_tracing() == false);
  bool was_on  = (simulation->projections_tracing() == true);
  bool turn_on = simulation->
    projections_schedule_on()->write_this_cycle(cycle_,time_);
  bool turn_off = simulation->
    projections_schedule_off()->write_this_cycle(cycle_,time_);

  if (was_off && turn_on) {

    simulation->monitor()->print
      ("Performance","turning projections logging ON\n");

    simulation->set_projections_tracing(true);

    traceBegin();

  } else if (was_on && turn_off) {

    simulation->monitor()->print
      ("Performance","turning projections logging OFF\n");

    simulation->set_projections_tracing(false);

    traceEnd();

  }
#endif

  stopping_balance_();

}

//----------------------------------------------------------------------

void CommBlock::stopping_balance_()
{

  int balance_interval =  simulation()->config()->balance_interval;
  
  if (balance_interval && ((cycle_ % balance_interval) == 0)) {
    VERBOSE("balance_enter");
    simulation()->set_phase (phase_balance);
    AtSync();
  } else {
    
    //    control_sync(phase_stopping_exit,"none",true,__FILE__,__LINE__);
    //    stopping_exit_();
    control_next();
  }
}

//----------------------------------------------------------------------

void CommBlock::ResumeFromSync()
{
  VERBOSE("balance_exit");
  //  control_sync(phase_stopping_exit,"none",true,__FILE__,__LINE__);
  control_next();
}

//----------------------------------------------------------------------

void CommBlock::exit_()
{
  if (index_.is_root()) {

    if (simulation())
      simulation()->performance_write();

    proxy_main.p_exit(1);
  }
}
