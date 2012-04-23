// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulationMpi user-dependent class member functions

#ifndef CONFIG_USE_CHARM

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationMpi::EnzoSimulationMpi
(
 const char * parameter_file,
 const GroupProcess * group_process ) throw ()
  : SimulationMpi(parameter_file,group_process)
{
}

//----------------------------------------------------------------------

EnzoSimulationMpi::~EnzoSimulationMpi() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::initialize() throw()
{
  SimulationMpi::initialize();
  EnzoBlock::initialize(parameters_,field_descr());
}

//----------------------------------------------------------------------

const Factory * EnzoSimulationMpi::factory() const throw()
{ 
  if (! factory_) factory_ = new EnzoFactory;
  return factory_;
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::run() throw()
{
  
  DEBUG("EnzoSimulationMpi::run()");

#ifdef CONFIG_USE_MPI
  ReduceMpi    reduce(group_process_);
#else
  ReduceSerial reduce(group_process_);
#endif

  Problem * problem = Simulation::problem();

  //--------------------------------------------------
  // INITIALIZE FIELDS
  //--------------------------------------------------

  DEBUG("EnzoSimulationMpi::run() Initial");

  Initial * initial = problem->initial();

  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      initial->enforce(hierarchy_,field_descr_,block);

      block->set_cycle(cycle_);
      block->set_time(time_);

      block->initialize();

    }
  }

  //--------------------------------------------------
  // REFRESH GHOST ZONES AND ENFORCE BOUNDARY CONDITIONS
  //--------------------------------------------------

  DEBUG("EnzoSimulationMpi::run() refresh, Boundary ");

  double lower[3];
  hierarchy_->lower(&lower[0], &lower[1], &lower[2]);
  double upper[3];
  hierarchy_->upper(&upper[0], &upper[1], &upper[2]);

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      bool is_boundary [3][2];

      block->is_on_boundary (lower,upper,is_boundary);

      refresh_ghost_(block,patch,is_boundary);

      update_boundary_(block,is_boundary);

    }
  }

  //--------------------------------------------------
  // PERFORM SCHEDULED OUTPUT
  //--------------------------------------------------


  DEBUG("EnzoSimulationMpi::run() Output");

  scheduled_output();
  
  //--------------------------------------------------
  // EVALUATE INITIAL STOPPING CRITERIA
  //--------------------------------------------------

  Stopping * stopping = problem->stopping();

  int stop_hierarchy = true;

  DEBUG("EnzoSimulationMpi::run() Stopping");

  while ((patch = ++it_patch)) {

    int    stop_patch  = true;

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      int cycle = block->cycle();
      int time  = block->time();

      int stop_block = stopping->complete(cycle,time);

      stop_patch = stop_patch && stop_block;
    }
    stop_hierarchy = stop_hierarchy && stop_patch;
  }

  //======================================================================
  // BEGIN MAIN LOOP
  //======================================================================

  while (! stop_hierarchy) {

    //--------------------------------------------------
    // Determine timestep
    //--------------------------------------------------


    Timestep * timestep = problem->timestep();
    Stopping * stopping = problem->stopping();

    double dt_hierarchy = std::numeric_limits<double>::max();

    // Accumulate Patch-local timesteps

    while ((patch = ++it_patch)) {

      double dt_patch = std::numeric_limits<double>::max();

      ItBlock it_block(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++it_block)) {

	double dt_block   = timestep->compute(field_descr_,block);
	double time_block = block->time();

	// Reduce timestep to coincide with scheduled output

	int index_output=0;
	while (Output * output = problem->output(index_output++)) {
	  dt_block = output->update_timestep(time_block,dt_block);
	}

	// Reduce timestep to coincide with end of simulation

	double time_stop = stopping->stop_time();
	dt_block = MIN (dt_block, (time_stop - time_block));

	dt_patch = MIN(dt_patch,dt_block);
      }
      dt_hierarchy = MIN(dt_hierarchy, dt_patch);
    }

    dt_hierarchy = reduce.reduce_double(dt_hierarchy,reduce_op_min);

    // Set Simulation timestep

    dt_ = dt_hierarchy;

    // Set Block timesteps

    while ((patch = ++it_patch)) {

      ItBlock it_block(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++it_block)) {
	block->set_dt    (dt_);
      }
    }

    ASSERT("EnzoSimulation::run", "dt == 0", dt_ != 0.0);

    monitor_output();

    //--------------------------------------------------
    // REFRESH GHOST ZONES AND ENFORCE BOUNDARY CONDITIONS
    //--------------------------------------------------

    while ((patch = ++it_patch)) {

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	bool is_boundary [3][2];

	block->is_on_boundary (lower,upper,is_boundary);

	refresh_ghost_(block,patch,is_boundary);

	update_boundary_(block,is_boundary);

      }
    }

    //--------------------------------------------------
    // APPLY NUMERICAL METHODS AND EVALUATE STOPPING CRITERIA
    //--------------------------------------------------

    stop_hierarchy = true;

    while ((patch = ++it_patch)) {

      int stop_patch = true;

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	int index_method = 0;
	while (Method * method = problem->method(index_method++)) {
	  method -> compute_block (field_descr_,block);
	}

	int   cycle_block = block->cycle();
	block->set_cycle (++cycle_block);

	double time_block = block->time();
	double   dt_block = block->dt();
	block->set_time  (time_block + dt_block);

	// Global cycle and time reduction
	
	int stop_block = stopping->complete(cycle_block,time_block);
	
	// Update stopping criteria for patch

	stop_patch = stop_patch && stop_block;

      } // (block = ++it_block)

      // Update stopping criteria for hierarchy

      stop_hierarchy = stop_hierarchy && stop_patch;

    } // (patch = ++it_patch)

    stop_hierarchy = reduce.reduce_int(stop_hierarchy,reduce_op_land);

    cycle_ ++;
    time_ += dt_hierarchy;

    // Perform any scheduled output

    scheduled_output();

  } // while (! stop_hierarchy)

  //======================================================================
  // END MAIN LOOP
  //======================================================================

  monitor_output();
  performance_output(performance_simulation_);

}

#endif /* ! CONFIG_USE_CHARM */
