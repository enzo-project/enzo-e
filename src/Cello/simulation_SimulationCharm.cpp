// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_SimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of SimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "simulation.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#include "simulation.hpp"


//----------------------------------------------------------------------

SimulationCharm::SimulationCharm
(
 const char         parameter_file[],
 int                n) throw ()
  : Simulation(parameter_file, n)
{
  // derived class should call initialize()
}

//----------------------------------------------------------------------

SimulationCharm::~SimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void SimulationCharm::initialize() throw()
{
  Simulation::initialize();
}

//----------------------------------------------------------------------

void SimulationCharm::run() throw()
{
  c_initial();
}

//----------------------------------------------------------------------

void SimulationCharm::s_initialize()
{

  if (patch_counter_.remaining() == 0) {
    DEBUG("Calling run()");
    run();
  }
  DEBUG("End s_initialize()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 //   ItPatch it_patch(hierarchy_);
//   Patch * patch;
//   while (( patch = ++it_patch )) {
//     CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
//     DEBUG1("proxy_patch = %p",proxy_patch);
//     DEBUG1("local patch = %p",proxy_patch->ckLocal());
//   }

}

//----------------------------------------------------------------------

void SimulationCharm::s_patch(CkCallback callback)
{
  DEBUG("s_patch");
  if (patch_counter_.remaining() == 0) {
    callback.send();
  }
}

//----------------------------------------------------------------------

void SimulationCharm::s_initial()
{
  DEBUG("s_initial");
  if (patch_counter_.remaining() == 0) {
    DEBUG ("SimulationCharm::s_initial() calling c_refresh()");
    c_refresh();
  } else  DEBUG ("SimulationCharm::s_initial() skipping c_refresh()");

}

//----------------------------------------------------------------------

void SimulationCharm::c_refresh()
{
  DEBUG ("SimulationCharm::c_refresh()");
  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while (( patch = ++it_patch )) {
    CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
    proxy_patch->p_refresh();
  }
}

//----------------------------------------------------------------------

void SimulationCharm::c_compute()
{
  //--------------------------------------------------
  // Stopping
  //--------------------------------------------------

  DEBUG1 ("SimulationCharm::c_compute() stop_ = %d",stop_);

  if (stop_) {
    DEBUG0;
    
    performance_output(performance_simulation_);

    DEBUG("Calling p_exit");
    proxy_main.p_exit(CkNumPes());
    DEBUG0;

  } else {
    DEBUG0;

    //--------------------------------------------------
    // Compute
    //--------------------------------------------------

    ItPatch it_patch(hierarchy_);
    Patch * patch;
    while (( patch = ++it_patch )) {
      CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
      DEBUG3("cycle %d time %f dt %f",
	     cycle_,time_,dt_);
      proxy_patch->p_compute(cycle_, time_, dt_);
    }
    DEBUG0;
  }
    DEBUG0;

}

//----------------------------------------------------------------------

void SimulationCharm::p_perf_output_min(CkReductionMsg * msg)
{
  // Collect minimum values

  double * reduce = (double * )msg->getData();
  for (size_t i=0; i<num_perf_; i++) {
    perf_min_[i] = reduce[i];
  }
  delete msg;

  // Then reduce maximum values

  CkCallback callback (CkIndex_SimulationCharm::p_perf_output_max(NULL),thisProxy);
  contribute( num_perf_*sizeof(double), perf_val_, 
	      CkReduction::max_double, callback);

}

//----------------------------------------------------------------------

void SimulationCharm::p_perf_output_max(CkReductionMsg * msg)
{
  // Collect maximum values

  double * reduce = (double * )msg->getData();
  for (size_t i=0; i<num_perf_; i++) {
    perf_max_[i] = reduce[i];
  }
  delete msg;

  // Finally reduce sum values

  CkCallback callback (CkIndex_SimulationCharm::p_perf_output_sum(NULL),thisProxy);
  contribute( num_perf_*sizeof(double), perf_val_, 
	      CkReduction::sum_double, callback);


}

//----------------------------------------------------------------------

void SimulationCharm::p_perf_output_sum(CkReductionMsg * msg)
{
  // Collect summed values

  double * reduce = (double * )msg->getData();
  for (size_t i=0; i<num_perf_; i++) {
    perf_sum_[i] = reduce[i];
  }
  delete msg;

  // Display performance output
  output_performance_();
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
