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
  : Simulation(parameter_file, n),
    block_sync_(0)
{
  TRACE("SimulationCharm::SimulationCharm");

}

//----------------------------------------------------------------------

SimulationCharm::~SimulationCharm() throw()
{
  TRACE("SimulationCharm::~SimulationCharm()");
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------

void SimulationCharm::performance_output()
{
  //  Simulation::performance_output();

  int num_regions  = performance_->num_regions();
  int num_counters =  performance_->num_counters();

  int n = num_regions * num_counters + 1;

  long long * counters_long_long = new long long [n];
  long *      counters_long = new long [n];

  for (int ir = 0; ir < num_regions; ir++) {
    performance_->region_counters(ir,counters_long_long);
    for (int ic = 0; ic < num_counters; ic++) {
      int index_counter = ir+num_regions*ic;
      counters_long[index_counter] = 
	(long) counters_long_long[index_counter];
    }
  }

  counters_long[n-1] = block_sync_.stop(); // number of CommBlocks

  CkCallback callback (CkIndex_SimulationCharm::p_performance_reduce(NULL), 
		       thisProxy);
  contribute (n*sizeof(long),
	      counters_long,CkReduction::sum_long,callback);
  delete [] counters_long;
  delete [] counters_long_long;

}

//----------------------------------------------------------------------

void SimulationCharm::p_performance_reduce(CkReductionMsg * msg)
{
  int num_regions  = performance_->num_regions();
  int num_counters =  performance_->num_counters();

  int n = num_regions * num_counters + 1;

  long long * counters_long_long = new long long [n];
  long *      counters_long = (long * )msg->getData();

  for (int ir = 0; ir < num_regions; ir++) {
    performance_->region_counters(ir,counters_long_long);
    for (int ic = 0; ic < num_counters; ic++) {
      int perf_counter = performance_->index_to_id(ic);
      int index_counter = ir+num_regions*ic;
      monitor_->print("Performance","%s %s %ld",
		      performance_->region_name(ir).c_str(),
		      performance_->counter_name(perf_counter).c_str(),
		      counters_long[ic]);  
    }
  }

  monitor_->print("Performance","simulation num-blocks %d",
		  counters_long[n-1]);

  delete [] counters_long_long;
  delete msg;

}
// void SimulationCharm::run() throw()
// {
//   TRACE("SimulationCharm::run()");
//   initial();
// }

//----------------------------------------------------------------------
//======================================================================

#endif /* CONFIG_USE_CHARM */
