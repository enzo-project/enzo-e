// See LICENSE_CELLO file for license and copyright information

/// @file     control_stopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
/// @ingroup  Control
///
///    STOPPING
///
///    Block::stopping()
///       update_boundary_()
///       compute dt
///       compute stopping
///       contribute( >>>>> Block::r_output() >>>>> )

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_STOPPING

#ifdef DEBUG_STOPPING
#   define TRACE_STOPPING(A)					\
  CkPrintf ("%d %s:%d %s TRACE %s\n",					\
	    CkMyPe(),__FILE__,__LINE__,name_.c_str(),A);
#else
#   define TRACE_STOPPING(A) ;
#endif


//----------------------------------------------------------------------

void Block::stopping_enter_()
{
  stopping_begin_();
}

//----------------------------------------------------------------------

void Block::stopping_begin_()
{

  TRACE_STOPPING("Block::stopping_begin_");

  Simulation * simulation = cello::simulation();

  simulation->set_phase(phase_stopping);

  int stopping_interval = cello::config()->stopping_interval;

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

    // Reduce to find Block array minimum dt and stopping criteria

    double min_reduce[2];

    min_reduce[0] = dt_block;
    min_reduce[1] = stop_block ? 1.0 : 0.0;

    CkCallback callback (CkIndex_Block::r_stopping_compute_timestep(NULL),
			 thisProxy);

#ifdef TRACE_CONTRIBUTE    
    CkPrintf ("%s %s:%d DEBUG_CONTRIBUTE\n",
	      name().c_str(),__FILE__,__LINE__); fflush(stdout);
#endif    
    contribute(2*sizeof(double), min_reduce, CkReduction::min_double, callback);

  } else {

    stopping_balance_();

  }

}

//----------------------------------------------------------------------

void Block::r_stopping_compute_timestep(CkReductionMsg * msg)
{
  performance_start_(perf_stopping);
  
  TRACE_STOPPING("Block::r_stopping_compute_timestep");
  
  ++age_;

  double * min_reduce = (double * )msg->getData();

  dt_   = min_reduce[0];
  stop_ = min_reduce[1] == 1.0 ? true : false;

  delete msg;

  Simulation * simulation = cello::simulation();

  dt_ *= Method::courant_global;
  
  set_dt   (dt_);
  set_stop (stop_);

  simulation->set_dt(dt_);
  simulation->set_stop(stop_);

#ifdef CONFIG_USE_PROJECTIONS
  bool was_off = (simulation->projections_tracing() == false);
  bool was_on  = (simulation->projections_tracing() == true);
  Schedule * schedule_on = simulation->projections_schedule_on();
  Schedule * schedule_off = simulation->projections_schedule_off();
  bool turn_on  = schedule_on  ? schedule_on->write_this_cycle(cycle_,time_) : false;
  bool turn_off = schedule_off ? schedule_off->write_this_cycle(cycle_,time_) : false;

  static bool active = false;
  if (!active && turn_on) {
    active = true;
    simulation->monitor()->print
      ("Performance","turning projections logging ON\n");

    simulation->set_projections_tracing(true);

    traceBegin();

  } else if (active && turn_off) {
    active = false;

    simulation->monitor()->print
      ("Performance","turning projections logging OFF\n");

    simulation->set_projections_tracing(false);

    traceEnd();

  }
#endif

  stopping_balance_();

  performance_stop_(perf_stopping);
}

//----------------------------------------------------------------------

void Block::stopping_balance_()
{
  TRACE_STOPPING("Block::stopping_balance_");

  Schedule * schedule = cello::simulation()->schedule_balance();

  bool do_balance = (schedule && 
		     schedule->write_this_cycle(cycle_,time_));

  if (do_balance) {

    const std::string balance_type = cello::config()->balance_type;

    if (balance_type == "cello") {

      // See EnzoMethodBalance

    } else if (balance_type == "charm") {

      // Charm++-controlled load balancing
      if (index_.is_root())
        cello::monitor()->print ("Balance","starting load balance step");
    }

    CkCallback callback = CkCallback
      (CkIndex_Block::r_stopping_load_balance(nullptr),
       proxy_array());
    adapt_ready_ = true;
    contribute(callback);
  } else {

    stopping_exit_();

  }
}

//----------------------------------------------------------------------

void Block::stopping_load_balance_()
{
  performance_start_(perf_stopping);
  TRACE_STOPPING("load_balance begin");
  cello::simulation()->set_phase (phase_balance);

  // Monitor * monitor = simulation()->monitor();
  // int mode_saved = monitor->mode();
  // monitor->set_mode(monitor_mode_all);
  // if (index().is_root()) monitor->print ("Balance","BEGIN");
  // monitor->set_mode(mode_saved);

  AtSync();
  performance_stop_(perf_stopping);
}

//----------------------------------------------------------------------

void Block::ResumeFromSync()
{
  // Monitor * monitor = simulation()->monitor();
  // int mode_saved = monitor->mode();
  // monitor->set_mode(monitor_mode_all);
  // if (index().is_root()) monitor->print ("Balance","END");
  // monitor->set_mode(mode_saved);

  TRACE_STOPPING("load_balance exit");

  stopping_exit_();

}

//----------------------------------------------------------------------

void Block::exit_()
{

  TRACE_STOPPING("Block::exit_");
  const int in = cello::index_static();
  if (index().is_root()) {
    if (DataMsg::counter[in] != 0) {
      CkPrintf ("%d Block::exit_() DataMsg::counter = %ld != 0\n",
		CkMyPe(),DataMsg::counter[in]);
      CkPrintf ("%d Block::exit_() ParticleData::counter = %ld != 0\n",
		CkMyPe(),ParticleData::counter[in]);
    }
    if (FieldFace::counter[in] != 0) {
      CkPrintf ("%d Block::exit_() FieldFace::counter = %ld != 0\n",
		CkMyPe(),FieldFace::counter[in]);
    }
    if (MsgCoarsen::counter[in] != 0) {
      CkPrintf ("%d Block::exit_() MsgCoarsen::counter = %ld != 0\n",
		CkMyPe(),MsgCoarsen::counter[in]);
    }
    if (MsgInitial::counter[in] != 0) {
      CkPrintf ("%d Block::exit_() MsgInitial::counter = %ld != 0\n",
		CkMyPe(),MsgInitial::counter[in]);
    }
    if (MsgOutput::counter[in] != 0) {
      CkPrintf ("%d Block::exit_() MsgOutput::counter = %ld != 0\n",
		CkMyPe(),MsgOutput::counter[in]);
    }
    if (MsgRefine::counter[in] != 0) {
      CkPrintf ("%d Block::exit_() MsgRefine::counter = %ld != 0\n",
		CkMyPe(),MsgRefine::counter[in]);
    }
    if (MsgRefresh::counter[in] != 0) {
      CkPrintf ("%d Block::exit_() MsgRefresh::counter = %ld != 0\n",
		CkMyPe(),MsgRefresh::counter[in]);
    }
  }
  if (index_.is_root()) {
    proxy_main.p_exit(1);
  }
}
