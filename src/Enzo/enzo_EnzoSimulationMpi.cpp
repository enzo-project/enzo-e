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

  // get hierarchy extents for later block boundary checks

  double lower[3];
  hierarchy_->lower(&lower[0], &lower[1], &lower[2]);
  double upper[3];
  hierarchy_->upper(&upper[0], &upper[1], &upper[2]);

  //--------------------------------------------------
  // INITIALIZE FIELDS
  //--------------------------------------------------

  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      initial_->enforce(hierarchy_,field_descr_,block);

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      enzo_block->set_cycle(cycle_);
      enzo_block->set_time(time_);

      enzo_block->initialize();

    }
  }

  while ((patch = ++it_patch)) {

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

      bool boundary [3][2];

      // is_block_on_boundary_ (block,boundary);
      block->is_on_boundary (lower,upper,boundary);

      update_boundary_(block,boundary);

      refresh_ghost_(block,patch,boundary);

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

	block->set_cycle(cycle_);
	block->set_time (time_);

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

    stop_hierarchy = true;

    //--------------------------------------------------
    // Refresh ghosts and boundary
    //--------------------------------------------------

    while ((patch = ++it_patch)) {

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	bool boundary [3][2];

	block->is_on_boundary (lower,upper,boundary);

	// Refresh ghost zones
	refresh_ghost_(block,patch,boundary);

	// Update boundary conditions
	update_boundary_(block,boundary);

      }
    }

    //--------------------------------------------------
    // Apply numerical methods
    //--------------------------------------------------

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

	  double lower[3];
	  double upper[3];
	  block->lower(lower+0,lower+1,lower+2);
	  block->upper(upper+0,upper+1,upper+2);
  char buffer[10];
  sprintf (buffer,"%03d-A",cycle_);
  block->field_block()->print(field_descr_,buffer,lower,upper);

	  method -> compute_block (field_descr_,block,time_block,dt_block);

  sprintf (buffer,"%03d-B",cycle_);
  block->field_block()->print(field_descr_,buffer,lower,upper);

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

//----------------------------------------------------------------------

void EnzoSimulationMpi::update_boundary_ (Block * block, bool boundary[3][2]) throw()
{
  // Update boundary conditions

  if (dimension_ >= 1) {
    if (boundary[axis_x][face_lower]) 
      boundary_->enforce(field_descr_,block,face_lower,axis_x);
    if (boundary[axis_x][face_upper]) 
      boundary_->enforce(field_descr_,block,face_upper,axis_x);
  }
  if (dimension_ >= 2) {
    if (boundary[axis_y][face_lower]) 
      boundary_->enforce(field_descr_,block,face_lower,axis_y);
    if (boundary[axis_y][face_upper]) 
      boundary_->enforce(field_descr_,block,face_upper,axis_y);
  }
  if (dimension_ >= 3) {
    if (boundary[axis_z][face_lower]) 
      boundary_->enforce(field_descr_,block,face_lower,axis_z);
    if (boundary[axis_z][face_upper]) 
      boundary_->enforce(field_descr_,block,face_upper,axis_z);
  }
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::refresh_ghost_ 
(
 Block * block, 
 Patch * patch, 
 bool    boundary[3][2]
 ) throw()
{
  // Refresh ghost zones

  int ibx,iby,ibz;
  block->index_patch(&ibx,&iby,&ibz);

  bool periodic = boundary_->is_periodic();

  if (dimension_ >= 1) {
    if (ibx % 2 == 0) {
      if (! boundary[axis_x][face_lower] || periodic) 
	block->refresh_ghosts(field_descr_,patch,-1,0,0);
      if (! boundary[axis_x][face_upper] || periodic) 
	block->refresh_ghosts(field_descr_,patch,+1,0,0);
    } else {
      if (! boundary[axis_x][face_upper] || periodic) 
	block->refresh_ghosts(field_descr_,patch,+1,0,0);
      if (! boundary[axis_x][face_lower] || periodic) 
	block->refresh_ghosts(field_descr_,patch,-1,0,0);
    }
  }
  if (dimension_ >= 2) {
    if (iby % 2 == 0) {
      if (! boundary[axis_y][face_lower] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,-1,0);
      if (! boundary[axis_y][face_upper] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,+1,0);
    } else {
      if (! boundary[axis_y][face_upper] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,+1,0);
      if (! boundary[axis_y][face_lower] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,-1,0);
    }

  }
  if (dimension_ >= 3) {
    if (ibz % 2 == 0) {
      if (! boundary[axis_z][face_lower] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,0,-1);
      if (! boundary[axis_z][face_upper] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,0,+1);
    } else {
      if (! boundary[axis_z][face_upper] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,0,+1);
      if (! boundary[axis_z][face_lower] || periodic) 
	block->refresh_ghosts(field_descr_,patch,0,0,-1);
    }
  }
}

#endif /* ! CONFIG_USE_CHARM */
