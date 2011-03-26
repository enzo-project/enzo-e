// $Id: enzo_EnzoSimulationCharm.cpp 2116 2011-03-17 22:11:13Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "enzo.hpp"

extern CProxy_Main mainProxy;

//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 const char parameter_file[],
 int        string_length,
 int        index) throw ()
  : EnzoSimulation(parameter_file, new GroupProcessCharm, index)
{

  //  Monitor::instance() -> set_active (CkMyPe() == 0);

  initialize();

  run();

  mainProxy.enzo_exit(index);

}

//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 const char * parameter_file,
 GroupProcess * group_process) throw ()
  : EnzoSimulation(parameter_file,group_process)
{
}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::run() throw()
{
  
  Monitor * monitor = Monitor::instance();

  Performance performance;

  performance.start();
  
  //--------------------------------------------------
  // INITIALIZE FIELDS
  //--------------------------------------------------

  ItPatch it_patch(mesh_);
  Patch * patch;

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      initial_->compute(block);

      boundary_->enforce(block);

      // Initialize Block attributes 
      // (REQUIRED HERE INSTEAD OF CONSTRUCTOR SINCE REQUIRES extents_[])

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      enzo_block->initialize();

    }
  }

  // Perform any scheduled output

  for (size_t i=0; i<output_list_.size(); i++) {
    Output * output = output_list_[i];
    output->scheduled_write(mesh_,cycle_,time_);
  }

  //--------------------------------------------------
  // INITIAL STOPPING CRITERIA TEST
  //--------------------------------------------------

  Reduce * reduce_mesh = mesh_->group()->create_reduce ();

  int    stop_patch  = true;

  while ((patch = ++it_patch)) {

    Reduce * reduce_patch = patch->group()->create_reduce ();
    int    stop_block  = true;

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      int cycle_block   = enzo_block->CycleNumber;
      double time_block =  enzo_block->Time;

      stop_block = stop_block && (stopping_->complete(cycle_block,time_block));
    }

    int stop_reduce = reduce_patch->reduce_int (stop_block, reduce_op_land);

    stop_patch  = stop_patch && stop_reduce;
		      
    delete reduce_patch;

  }

  int stop_mesh = reduce_mesh->reduce_int (stop_patch, reduce_op_land);
  
  //======================================================================
  // BEGIN MAIN LOOP
  //======================================================================

  while (! stop_mesh) {

    monitor->print("[Simulation %d] cycle %04d time %15.12f",
		   index_, cycle_,time_);

    //--------------------------------------------------
    // Determine timestep and dump output
    //--------------------------------------------------

    double dt_patch = std::numeric_limits<double>::max();

    // Accumulate Patch-local timesteps

    while ((patch = ++it_patch)) {

      Reduce * reduce_patch = patch->group()->create_reduce ();

      double dt_block = std::numeric_limits<double>::max();

      ItBlock it_block(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++it_block)) {

	// Enforce boundary conditions
	boundary_->enforce(block);

	// Accumulate Block-local dt

	dt_block = MIN(dt_block,timestep_->compute(block));

	// Reduce timestep to coincide with scheduled output if needed

	EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

	double time = enzo_block->Time;

	for (size_t i=0; i<output_list_.size(); i++) {
	  Output * output = output_list_[i];
	  dt_block = output->update_timestep(time,dt_block);
	}

      } // ( block = ++it_block )

      double dt_reduce = reduce_patch->reduce_double (dt_block, reduce_op_min);

      dt_patch = MIN(dt_patch, dt_reduce);

      delete reduce_patch;

    } // ( patch = ++it_patch )

    double dt_mesh  = reduce_mesh->reduce_double(dt_patch,reduce_op_min);

    ASSERT("EnzoSimulation::run", "dt == 0", dt_mesh != 0.0);

    //--------------------------------------------------
    // Apply the methods
    //--------------------------------------------------

    // Initialize reductions

    int  stop_patch = true;

    while ((patch = ++it_patch)) {

      Reduce * reduce_patch = patch->group()->create_reduce ();

      int stop_block  = true;

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	// Refresh ghosts

	block->refresh_ghosts();

	// Initialize enzo_block aliases (may me modified)

	EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

	int &    cycle_block        = enzo_block->CycleNumber;
	enzo_float & time_block     = enzo_block->Time;
	enzo_float & dt_block       = enzo_block->dt;
	enzo_float & old_time_block = enzo_block->OldTime;

	// UNIFORM TIMESTEP OVER ALL BLOCKS IN MESH

	double dt_stop = (stopping_->stop_time() - dt_block);

	dt_block = MIN(dt_mesh, dt_stop);

	// Loop through methods

	for (size_t i = 0; i < method_list_.size(); i++) {

	  Method * method = method_list_[i];

	  method -> compute_block (block,time_block,dt_block);

	} // method

	// Update enzo_block values

	old_time_block  = time_block;
	time_block     += dt_block;
	cycle_block    += 1;

	// Global cycle and time reduction
	
	stop_block = stop_block 
	  &&  (stopping_->complete(cycle_block,time_block));

      } // (block = ++it_block)

      int stop_reduce = reduce_patch->reduce_int (stop_block, reduce_op_land);

      stop_patch  = stop_patch && stop_reduce;

      delete reduce_patch;
    } // (patch = ++it_patch)

    stop_mesh  = reduce_mesh->reduce_int (stop_patch, reduce_op_land);

    cycle_ ++;
    time_ += dt_mesh;

    // Perform any scheduled output

    for (size_t i=0; i<output_list_.size(); i++) {
      Output * output = output_list_[i];
      output->scheduled_write(mesh_,cycle_,time_);
    }
  } // ! stop

  //======================================================================
  // END MAIN LOOP
  //======================================================================

  monitor->print("[Simulation %d] cycle %04d time %15.12f", 
		 index_,
		 cycle_,
		 time_);

  performance.stop();

  performance.print();


}

//======================================================================

#endif /* CONFIG_USE_CHARM */
