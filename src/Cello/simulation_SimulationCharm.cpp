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

  WARNING("SimulationCharm::initialize()",
	  "Calling StartLB for debugging load balancing()");

}

//----------------------------------------------------------------------

void SimulationCharm::run() throw()
{
  c_initial();
}

//----------------------------------------------------------------------

void SimulationCharm::s_initialize()
{
  if (patch_loop_.done()) run();
}

//----------------------------------------------------------------------

void SimulationCharm::s_patch(CkCallback callback)
{
  if (patch_loop_.done()) callback.send();
}

//----------------------------------------------------------------------

void SimulationCharm::c_refresh()
{
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
  if (stop_) {
    
    performance_output(performance_simulation_);

    proxy_main.p_exit(CkNumPes());

  } else {

    ItPatch it_patch(hierarchy_);
    Patch * patch;
    while (( patch = ++it_patch )) {
      CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
      proxy_patch->p_compute(cycle_, time_, dt_);
    }
  }
}

//----------------------------------------------------------------------

void SimulationCharm::p_performance_min(CkReductionMsg * msg)
{
  // Collect minimum values

  double * reduce = (double * )msg->getData();
  for (size_t i=0; i<num_perf_; i++) {
    perf_min_[i] = reduce[i];
  }
  delete msg;

  // Then reduce maximum values

  CkCallback callback (CkIndex_SimulationCharm::p_performance_max(NULL),thisProxy);
  contribute( num_perf_*sizeof(double), perf_val_, 
	      CkReduction::max_double, callback);

}

//----------------------------------------------------------------------

void SimulationCharm::p_performance_max(CkReductionMsg * msg)
{
  // Collect maximum values

  double * reduce = (double * )msg->getData();
  for (size_t i=0; i<num_perf_; i++) {
    perf_max_[i] = reduce[i];
  }
  delete msg;

  // Finally reduce sum values

  CkCallback callback (CkIndex_SimulationCharm::p_performance_sum(NULL),thisProxy);
  contribute( num_perf_*sizeof(double), perf_val_, 
	      CkReduction::sum_double, callback);


}

//----------------------------------------------------------------------

void SimulationCharm::p_performance_sum(CkReductionMsg * msg)
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
