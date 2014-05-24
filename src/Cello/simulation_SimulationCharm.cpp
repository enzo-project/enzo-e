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
    sync_output_begin_(0),
    sync_output_write_(0)
{
  TRACE("SimulationCharm::SimulationCharm");
}

//----------------------------------------------------------------------

SimulationCharm::~SimulationCharm() throw()
{
  TRACE("SimulationCharm::~SimulationCharm()");
#ifdef CELLO_DEBUG
  fclose (fp_debug_);
#endif
}

//----------------------------------------------------------------------

void SimulationCharm::insert_block() 
{
 
#ifdef CELLO_DEBUG
  PARALLEL_PRINTF ("%d: ++sync_output_begin_ %d %d\n",
		   CkMyPe(),sync_output_begin_.stop(),hierarchy()->num_blocks());
#endif
  hierarchy()->increment_block_count(1);
  ++sync_output_begin_;
  ++sync_output_write_;
}

//----------------------------------------------------------------------

void SimulationCharm::delete_block() 
{
#ifdef CELLO_DEBUG
  PARALLEL_PRINTF ("%d: --block_sync_ %d %d\n",
		   CkMyPe(),sync_output_begin_.stop(),
		   hierarchy()->num_blocks());
#endif
  hierarchy()->increment_block_count(-1);
  --sync_output_begin_;
  --sync_output_write_;
}

//----------------------------------------------------------------------

void SimulationCharm::p_monitor()
{
  monitor()-> print("", "-------------------------------------");
  monitor()-> print("Simulation", "cycle %04d", cycle_);
  monitor()-> print("Simulation", "time-sim %15.12f",time_);
  monitor()-> print("Simulation", "dt %15.12g", dt_);
  proxy_simulation.p_monitor_performance();
}

//----------------------------------------------------------------------

void SimulationCharm::monitor_performance()
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

  counters_long[n-1] = hierarchy()->num_blocks(); // number of CommBlocks

  // --------------------------------------------------
  // ENTRY: #1 SimulationCharm::monitor_performance() -> SimulationCharm::r_monitor_performance()
  // ENTRY: contribute()
  // --------------------------------------------------
  CkCallback callback (CkIndex_SimulationCharm::r_monitor_performance(NULL), 
		       thisProxy);
  contribute (n*sizeof(long), counters_long,CkReduction::sum_long,callback);
  // --------------------------------------------------

  delete [] counters_long;
  delete [] counters_long_long;

}

//----------------------------------------------------------------------

void SimulationCharm::r_monitor_performance(CkReductionMsg * msg)
{
  int nr  = performance_->num_regions();
  int nc =  performance_->num_counters();

  int n = nr * nc + 1;

  long *      counters_long = (long * )msg->getData();

  int index_region_cycle = performance_->region_index("cycle");

  for (int ir = 0; ir < nr; ir++) {
    for (int ic = 0; ic < nc; ic++) {
      int index_counter = ir+nr*ic;
      bool do_print = 
	(performance_->counter_type(ic) != counter_type_abs) ||
	(ir == index_region_cycle);
      if (do_print) {
	monitor()->print("Performance","%s %s %ld",
			performance_->region_name(ir).c_str(),
			performance_->counter_name(ic).c_str(),
			counters_long[index_counter]);
      }
    }
  }

  monitor()->print("Performance","simulation num-blocks %d",
		  counters_long[n-1]);

  Memory::instance()->reset_high();

  delete msg;

}

//----------------------------------------------------------------------

