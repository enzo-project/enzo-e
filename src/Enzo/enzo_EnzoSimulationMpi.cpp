// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @todo     Move all block initialization code into block constructor
/// @todo Move check for whether a block face is along a boundary
///           within Boundary::enforce()
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulationMpi user-dependent class member functions

#ifndef CONFIG_USE_CHARM

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationMpi::EnzoSimulationMpi
(
 const char * parameter_file,
 GroupProcess * group_process,
 int index ) throw ()
  : EnzoSimulation(parameter_file,group_process, index)
{
}

//----------------------------------------------------------------------

EnzoSimulationMpi::~EnzoSimulationMpi() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::run() throw()
{
  
  Performance performance;
#ifdef CONFIG_USE_MPI
  ReduceMpi    reduce(group_process_);
#else
  ReduceSerial reduce(group_process_);
#endif

  performance.start();
  
  //--------------------------------------------------
  // INITIALIZE FIELDS
  //--------------------------------------------------

  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      initial_->compute(field_descr_,block);

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      enzo_block->initialize(cycle_, time_);

    }
  }


  double lower_hierarchy[3];
  hierarchy_->lower(&lower_hierarchy[0],
		    &lower_hierarchy[1],
		    &lower_hierarchy[2]);

  double upper_hierarchy[3];
  hierarchy_->upper(&upper_hierarchy[0],
		    &upper_hierarchy[1],
		    &upper_hierarchy[2]);

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

	double lower_block[3];
	block->lower(&lower_block[axis_x],
		     &lower_block[axis_y],
		     &lower_block[axis_z]);
	double upper_block[3];
	block->upper(&upper_block[axis_x],
		     &upper_block[axis_y],
		     &upper_block[axis_z]);

	int ix,iy,iz;
	block->index_patch(&ix,&iy,&iz);
	bool lower_x = 
	  (cello::err_rel(lower_block[axis_x],lower_hierarchy[axis_x]) < 1e-7);
	bool lower_y = 
	  (cello::err_rel(lower_block[axis_y],lower_hierarchy[axis_y]) < 1e-7);
	bool lower_z = 
	  (cello::err_rel(lower_block[axis_z],lower_hierarchy[axis_z]) < 1e-7);
	bool upper_x = 
	  (cello::err_rel(upper_block[axis_x],upper_hierarchy[axis_x]) < 1e-7);
	bool upper_y = 
	  (cello::err_rel(upper_block[axis_y],upper_hierarchy[axis_y]) < 1e-7);
	bool upper_z = 
	  (cello::err_rel(upper_block[axis_z],upper_hierarchy[axis_z]) < 1e-7);

	if (dimension_ >= 1) {
	  if (lower_x) boundary_->enforce(field_descr_,block,face_lower,axis_x);
	  if (upper_x) boundary_->enforce(field_descr_,block,face_upper,axis_x);
	}
	if (dimension_ >= 2) {
	  if (lower_y) boundary_->enforce(field_descr_,block,face_lower,axis_y);
	  if (upper_y) boundary_->enforce(field_descr_,block,face_upper,axis_y);
	}
	if (dimension_ >= 3) {
	  if (lower_z) boundary_->enforce(field_descr_,block,face_lower,axis_z);
	  if (upper_z) boundary_->enforce(field_descr_,block,face_upper,axis_z);
	}

	// Refresh ghost zones

	if (dimension_ >= 1) {
	  if (! lower_x) block->refresh_ghosts(field_descr_,patch,face_lower,axis_x);
	  if (! upper_x) block->refresh_ghosts(field_descr_,patch,face_upper,axis_x);
	}
	if (dimension_ >= 2) {
	  if (! lower_y) block->refresh_ghosts(field_descr_,patch,face_lower,axis_y);
	  if (! upper_y) block->refresh_ghosts(field_descr_,patch,face_upper,axis_y);
	}
	if (dimension_ >= 3) {
	  if (! lower_z) block->refresh_ghosts(field_descr_,patch,face_lower,axis_z);
	  if (! upper_z) block->refresh_ghosts(field_descr_,patch,face_upper,axis_z);
	}
    }
  }
  // Perform any scheduled output

  scheduled_output();
  
  //--------------------------------------------------
  // INITIAL STOPPING CRITERIA TEST
  //--------------------------------------------------

  int stop_hierarchy = true;

  while ((patch = ++it_patch)) {

    int    stop_patch  = true;

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      int & cycle_block = enzo_block->CycleNumber;
      double time_block =  enzo_block->Time();

      int stop_block = stopping_->complete(cycle_block,time_block);

      stop_patch = stop_patch && stop_block;

    }

    stop_hierarchy = stop_hierarchy && stop_patch;

  }

  int ip = group_process_->rank();
   int np = group_process_->size();

  //======================================================================
  // BEGIN MAIN LOOP
  //======================================================================

  while (! stop_hierarchy) {

    //--------------------------------------------------
    // Determine timestep
    //--------------------------------------------------

    double dt_hierarchy = std::numeric_limits<double>::max();

    // Accumulate Patch-local timesteps

    while ((patch = ++it_patch)) {

      double dt_patch = std::numeric_limits<double>::max();

      ItBlock it_block(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++it_block)) {

	// Accumulate Block-local dt

	block->initialize(cycle_,time_);

	double dt_block = timestep_->compute(field_descr_,block);

	// Reduce timestep to coincide with scheduled output if needed

	double time_block = static_cast <EnzoBlock*> (block)->Time();

	for (size_t i=0; i<output_list_.size(); i++) {
	  dt_block = output_list_[i]->update_timestep(time_block,dt_block);
	}

	// Reduce timestep to coincide with end of simulation if needed

	dt_block = MIN (dt_block, (stopping_->stop_time() - time_block));

	// Update patch-level timestep

	dt_patch = MIN(dt_patch,dt_block);

      } // ( block = ++it_block )

      // Update hierarchy-level timestep

      dt_hierarchy = MIN(dt_hierarchy, dt_patch);

    } // ( patch = ++it_patch )

    dt_hierarchy = reduce.reduce_double(dt_hierarchy,reduce_op_min);

    dt_ = dt_hierarchy;

    ASSERT("EnzoSimulation::run", "dt == 0", dt_hierarchy != 0.0);

    monitor_output();

    //--------------------------------------------------
    // Apply the methods
    //--------------------------------------------------

    stop_hierarchy = true;


    while ((patch = ++it_patch)) {

      int stop_patch = true;

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	// Enforce boundary conditions ; refresh ghost zones

	double lower_block[3];
	block->lower(&lower_block[axis_x],
		     &lower_block[axis_y],
		     &lower_block[axis_z]);
	double upper_block[3];
	block->upper(&upper_block[axis_x],
		     &upper_block[axis_y],
		     &upper_block[axis_z]);

	int ix,iy,iz;
	block->index_patch(&ix,&iy,&iz);
	bool lower_x = 
	  (cello::err_rel(lower_block[axis_x],lower_hierarchy[axis_x]) < 1e-7);
	bool lower_y = 
	  (cello::err_rel(lower_block[axis_y],lower_hierarchy[axis_y]) < 1e-7);
	bool lower_z = 
	  (cello::err_rel(lower_block[axis_z],lower_hierarchy[axis_z]) < 1e-7);
	bool upper_x = 
	  (cello::err_rel(upper_block[axis_x],upper_hierarchy[axis_x]) < 1e-7);
	bool upper_y = 
	  (cello::err_rel(upper_block[axis_y],upper_hierarchy[axis_y]) < 1e-7);
	bool upper_z = 
	  (cello::err_rel(upper_block[axis_z],upper_hierarchy[axis_z]) < 1e-7);

	// Update boundary conditions

	if (dimension_ >= 1) {
	  if (lower_x) boundary_->enforce(field_descr_,block,face_lower,axis_x);
	  if (upper_x) boundary_->enforce(field_descr_,block,face_upper,axis_x);
	}
	if (dimension_ >= 2) {
	  if (lower_y) boundary_->enforce(field_descr_,block,face_lower,axis_y);
	  if (upper_y) boundary_->enforce(field_descr_,block,face_upper,axis_y);
	}
	if (dimension_ >= 3) {
	  if (lower_z) boundary_->enforce(field_descr_,block,face_lower,axis_z);
	  if (upper_z) boundary_->enforce(field_descr_,block,face_upper,axis_z);
	}

	// Refresh ghost zones

	if (dimension_ >= 1) {
	  if (! lower_x) block->refresh_ghosts(field_descr_,patch,face_lower,axis_x);
	  if (! upper_x) block->refresh_ghosts(field_descr_,patch,face_upper,axis_x);
	}
	if (dimension_ >= 2) {
	  if (! lower_y) block->refresh_ghosts(field_descr_,patch,face_lower,axis_y);
	  if (! upper_y) block->refresh_ghosts(field_descr_,patch,face_upper,axis_y);
	}
	if (dimension_ >= 3) {
	  if (! lower_z) block->refresh_ghosts(field_descr_,patch,face_lower,axis_z);
	  if (! upper_z) block->refresh_ghosts(field_descr_,patch,face_upper,axis_z);
	}


      }
    }

    while ((patch = ++it_patch)) {

      int stop_patch = true;

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	double lower_block[3];
	block->lower(&lower_block[axis_x],
		     &lower_block[axis_y],
		     &lower_block[axis_z]);
	double upper_block[3];
	block->upper(&upper_block[axis_x],
		     &upper_block[axis_y],
		     &upper_block[axis_z]);
	FieldBlock * field_block = block -> field_block();

	double * values = (double * ) 
	       field_block->field_unknowns(field_descr_, 0);

	EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

	int        & cycle_block    = enzo_block->CycleNumber;
	enzo_float   time_block     = enzo_block->Time();
	enzo_float & dt_block       = enzo_block->dt;
	enzo_float & old_time_block = enzo_block->OldTime;

	// UNIFORM TIMESTEP OVER ALL BLOCKS IN HIERARCHY

	dt_block = dt_hierarchy;

	// Loop through methods

	for (size_t i = 0; i < method_list_.size(); i++) {

	  Method * method = method_list_[i];

	  method -> compute_block (block,time_block,dt_block);

	}

	// Update enzo_block values

	old_time_block  = time_block;
	time_block     += dt_block;
	cycle_block    += 1;

	// Global cycle and time reduction
	
	int stop_block = stopping_->complete(cycle_block,time_block);
	
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

  performance.stop();

  performance.print(monitor_);

}

#endif /* ! CONFIG_USE_CHARM */
