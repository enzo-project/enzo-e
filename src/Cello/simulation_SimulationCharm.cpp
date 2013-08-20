// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_SimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of SimulationCharm user-dependent class member functions

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

  int nr  = performance_->num_regions();
  int nc =  performance_->num_counters();

  int n = nr * nc + 1;

  long long * counters_long_long = new long long [nc];
  long *      counters_long = new long [n];

  for (int ir = 0; ir < nr; ir++) {
    performance_->region_counters(ir,counters_long_long);
    for (int ic = 0; ic < nc; ic++) {
      int index_counter = ir+nr*ic;
      counters_long[index_counter] = 
	(long) counters_long_long[ic];
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
  int nr  = performance_->num_regions();
  int nc =  performance_->num_counters();

  int n = nr * nc + 1;

  long *      counters_long = (long * )msg->getData();

  for (int ir = 0; ir < nr; ir++) {
    for (int ic = 0; ic < nc; ic++) {
      int index_counter = ir+nr*ic;
      bool do_print = 
	(performance_->counter_type(ic) != counter_type_abs) ||
	(ir == 0);
	
      if (do_print) {
	monitor_->print("Performance","%s %s %ld",
			performance_->region_name(ir).c_str(),
			performance_->counter_name(ic).c_str(),
			counters_long[index_counter]);
      }
    }
  }

  monitor_->print("Performance","simulation num-blocks %d",
		  counters_long[n-1]);

  delete msg;

}

//----------------------------------------------------------------------

