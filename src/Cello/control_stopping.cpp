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

// #define TRACE_DT

// #define DEBUG_STATE

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
    ((state_->cycle() % stopping_interval) == 0) : false;

  if (stopping_reduce || state_->dt()==0.0) {

    // Compute dt_ik for block i and method k

    Problem * problem = simulation->problem();

    //    allocate reduction vector for stopping criteria plus method dt
    const int n = 1 + problem->num_methods();
    std::vector<double> min_reduce(n);

    // Evaluate local stopping criteria

    Stopping * stopping = problem->stopping();
    const int stop_block = stopping->complete(state_->cycle(),state_->time());
    min_reduce[0] = stop_block ? 1.0 : 0.0;

    // Evaluate dt for each method k
    for (int k=0; k<problem->num_methods(); k++) {
      min_reduce[k+1] = problem->method(k)->timestep(this);
    }

    CkCallback callback (CkIndex_Block::r_stopping_compute_timestep(NULL),
			 thisProxy);

#ifdef TRACE_CONTRIBUTE
    CkPrintf ("%s %s:%d DEBUG_CONTRIBUTE\n",
	      name().c_str(),__FILE__,__LINE__); fflush(stdout);
#endif
    contribute(n*sizeof(double), min_reduce.data(), CkReduction::min_double, callback);

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

  state_->set_stopping(min_reduce[0] == 1.0);

  // Compute timestep
  double dt_global = stopping_compute_global_dt_(min_reduce);

  stopping_update_method_state_(min_reduce,dt_global);

  delete msg;

  // Update Block state timestep
  state_->set_dt(dt_global);

  // Update simulation state to block state
  Simulation * simulation = cello::simulation();
  simulation->state()->set_dt      (state_->dt());
  simulation->state()->set_stopping(state_->stopping());

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

double Block::stopping_compute_global_dt_ (double min_reduce[])
{
  Simulation * simulation = cello::simulation();
  Problem * problem = simulation->problem();
  double dt_global = std::numeric_limits<double>::max();

  // compute minimum timestep dt_global over all methods
  for (int k=0; k<problem->num_methods(); k++) {
    const double dt_method = min_reduce[k+1];
    dt_global = std::min(dt_global,dt_method);
  }

  // Adjust timestep dt for global courant condition
  dt_global *= Method::courant_global;

  // adjust timestep dt to align with any scheduled output times
  double time_curr = state_->time();
  int index_output=0;
  while (Output * output = problem->output(index_output++)) {
    Schedule * schedule = output->schedule();
    dt_global = schedule->update_timestep(time_curr,dt_global);
  }

  // Reduce timestep to not overshoot final time from stopping criteria

  Stopping * stopping = problem->stopping();
  double time_stop = stopping->stop_time();

  dt_global = std::min (dt_global, (time_stop - time_curr));

  return dt_global;
}

//----------------------------------------------------------------------

void Block::stopping_update_method_state_(double min_reduce[], double dt_global)
{
  // update Method states for supercycling
  Simulation * simulation = cello::simulation();
  Problem * problem = simulation->problem();
#ifdef DEBUG_STATE
  // Write current state
  if (index().is_root()) {
    state()->print("update_method_state");
  }
#endif
  for (int k=0; k<problem->num_methods(); k++) {
    const double dt_method = min_reduce[k+1];
    const int max_super = problem->method(k)->max_supercycle();
    const double max_dt_method = dt_global*max_super;
    const double ratio = dt_method / dt_global;
    const int desired_super = int(std::floor(ratio));
    const int allowed_super = std::min(desired_super,max_super);
    State::MethodState & method_state = state_->method(k);
    const int step = method_state.step();
    const int num_steps = method_state.num_steps();
    bool update_state = (step >= num_steps);
    if (update_state) {
      method_state.set_time(state_->time());
      method_state.set_dt(allowed_super * dt_global);
      method_state.set_num_steps(allowed_super);
      method_state.set_step(0);
    }
  }
}

//----------------------------------------------------------------------

void Block::stopping_balance_()
{
  TRACE_STOPPING("Block::stopping_balance_");

  Schedule * schedule = cello::simulation()->schedule_balance();

  bool do_balance = (schedule && 
		     schedule->write_this_cycle(state_->cycle(),state_->time()));

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
