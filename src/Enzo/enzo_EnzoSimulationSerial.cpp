// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationSerial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move boundary conditions to Boundary object; remove bc_reflecting
/// @todo     Move output to Output object
/// @todo     Move timestep reductions into Timestep object
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulationSerial user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationSerial::EnzoSimulationSerial
(
 Parameters * parameters,
 GroupProcess * group_process) throw ()
  : Simulation(parameters,group_process)
{
}

//----------------------------------------------------------------------

EnzoSimulationSerial::~EnzoSimulationSerial() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

void EnzoSimulationSerial::initialize() throw()
{
  // Call initialize for Simulation base class
  Simulation::initialize();

  // Initialize enzo namespace variables
  enzo::initialize(parameters_);

}

//----------------------------------------------------------------------

void EnzoSimulationSerial::finalize() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationSerial::run() throw()
{
  
  Monitor * monitor = Monitor::instance();

  Performance performance;

  performance.start();
  
  //--------------------------------------------------
  // INITIALIZE FIELDS
  //--------------------------------------------------

  ItPatch itPatch(mesh_);
  Patch * patch;

  while ((patch = ++itPatch)) {

    ItBlock itBlock(patch);
    Block * block;

    while ((block = ++itBlock)) {

      initial_->compute(block);

      // Initialize Block attributes 
      // (REQUIRED HERE INSTEAD OF CONSTRUCTOR SINCE REQUIRES extents_[])

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      enzo_block->initialize();

      // Initial output (repeated below: remove when Output implement)

      block->field_block()->enforce_boundary(boundary_reflecting);
      output_images_(block, "enzo-p-%06d.%d.png",0,1);
    }
  }

  //--------------------------------------------------
  // INITIAL STOPPING CRITERIA TEST
  //--------------------------------------------------

  Reduce * reduce_mesh = mesh_->group()->create_reduce ();

  int    cycle_patch = std::numeric_limits<int>::max();
  int    stop_patch  = true;
  double time_patch  = std::numeric_limits<double>::max();

  while ((patch = ++itPatch)) {

    Reduce * reduce_patch = patch->group()->create_reduce ();

    int    cycle_block = std::numeric_limits<int>::max();
    int    stop_block  = true;
    double time_block  = std::numeric_limits<double>::max();

    ItBlock itBlock(patch);
    Block * block;

    while ((block = ++itBlock)) {

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      int cycle   = enzo_block->CycleNumber;
      double time = enzo_block->Time;

      cycle_block = MIN(cycle_block, cycle);
      time_block  = MIN(time_block,  time);
      stop_block  = stop_block && (stopping_->complete(cycle,time));
    }

    int   cycle_reduce = reduce_patch->reduce_int (cycle_block, reduce_op_min);
    double time_reduce = reduce_patch->reduce_double (time_block, reduce_op_min);
    int    stop_reduce = reduce_patch->reduce_int (stop_block, reduce_op_land);

    cycle_patch = MIN(cycle_patch, cycle_reduce);
    time_patch  = MIN(time_patch, time_reduce);
    stop_patch  = stop_patch && stop_reduce;
		      
    delete reduce_patch;

  }

  int   cycle_mesh = reduce_mesh->reduce_int (cycle_patch, reduce_op_min);
  double time_mesh = reduce_mesh->reduce_double (time_patch, reduce_op_min);
  int    stop_mesh = reduce_mesh->reduce_int (stop_patch, reduce_op_land);
  
  //======================================================================
  // BEGIN MAIN LOOP
  //======================================================================

  while (! stop_mesh) {

    monitor->print("cycle %04d time %15.12f", cycle_mesh,time_mesh);

    //--------------------------------------------------
    // Determine timestep and dump output
    //--------------------------------------------------

    double dt_patch = std::numeric_limits<double>::max();

    // Accumulate Patch-local timesteps

    while ((patch = ++itPatch)) {

      Reduce * reduce_patch = patch->group()->create_reduce ();

      double dt_block = std::numeric_limits<double>::max();

      ItBlock itBlock(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++itBlock)) {

	// Enforce boundary conditions
	block->field_block()->enforce_boundary(boundary_reflecting);

	// Output while we're here
	output_images_(block, "enzo-p-%06d.%d.png",cycle_mesh,0);

	// Accumulate Block-local dt

	dt_block = MIN(dt_block,timestep_->compute(block));

      } // ( block = ++itBlock )

      double dt_reduce = reduce_patch->reduce_double (dt_block, reduce_op_min);
      dt_patch = MIN(dt_patch, dt_reduce);

      delete reduce_patch;

    } // ( patch = ++itPatch )

    double dt_mesh  = reduce_mesh->reduce_double(dt_patch,reduce_op_min);

    ASSERT("EnzoSimulationSerial::run", "dt == 0", dt_mesh != 0.0);

    //--------------------------------------------------
    // Apply the methods
    //--------------------------------------------------

    // Initialize reductions

    int    cycle_patch = std::numeric_limits<int>::max();
    int    stop_patch  = true;
    double time_patch  = std::numeric_limits<double>::max();

    while ((patch = ++itPatch)) {

      Reduce * reduce_patch = patch->group()->create_reduce ();

      int    cycle_block = std::numeric_limits<int>::max();
      int    stop_block  = true;
      enzo_float time_block  = std::numeric_limits<enzo_float>::max();

      ItBlock itBlock(patch);
      Block * block;

      while ((block = ++itBlock)) {

	EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

	// Initialize enzo_block aliases (may me modified)

	int &    cycle        = enzo_block->CycleNumber;
	enzo_float & time     = enzo_block->Time;
	enzo_float & dt       = enzo_block->dt;
	enzo_float & old_time = enzo_block->OldTime;

	// UNIFORM TIMESTEP OVER ALL BLOCKS IN MESH

	double dt_stop = (stopping_->stop_time() - dt);
	dt = MIN(dt_mesh, dt_stop);

	// Loop through methods

	for (size_t i = 0; i < method_list_.size(); i++) {

	  Method * method = method_list_[i];

	  method -> compute_block (block,time,dt);

	} // method

	// Update enzo_block values

	cycle    += 1;
	old_time  = time;
	time     += dt;

	// Global cycle and time reduction
	
	cycle_block = MIN(cycle_block, cycle);
	time_block  = MIN(time_block,  time);
	stop_block  = stop_block && (stopping_->complete(cycle,time));

      } // (block = ++itBlock)

      int cycle_reduce = reduce_patch->reduce_int (cycle_block, reduce_op_min);
      double time_reduce = reduce_patch->reduce_double (time_block, reduce_op_min);
      int    stop_reduce = reduce_patch->reduce_int (stop_block, reduce_op_land);

      cycle_patch = MIN(cycle_patch, cycle_reduce);
      time_patch  = MIN(time_patch, time_reduce);
      stop_patch  = stop_patch && stop_reduce;

      delete reduce_patch;
    } // (patch = ++itPatch)

    cycle_mesh = reduce_mesh->reduce_int (cycle_patch, reduce_op_min);
    time_mesh  = reduce_mesh->reduce_double (time_patch, reduce_op_min);
    stop_mesh  = reduce_mesh->reduce_int (stop_patch, reduce_op_land);

  } // ! stop

  //======================================================================
  // END MAIN LOOP
  //======================================================================

  monitor->print("cycle %04d time %15.12f", cycle_mesh,time_mesh);

  //--------------------------------------------------
  // Final output dump
  //--------------------------------------------------

  while ((patch = ++itPatch)) {

    ItBlock itBlock(patch);
    Block * block;

    while ((block = ++itBlock)) {

      output_images_(block, "enzo-p-%06d.%d.png",cycle_mesh,1);

    }
  }

  performance.stop();

  performance.print();

}

//----------------------------------------------------------------------

void EnzoSimulationSerial::read() throw()
{
  INCOMPLETE("EnzoSimulationSerial::read");
}

//----------------------------------------------------------------------

void EnzoSimulationSerial::write() throw()
{
  INCOMPLETE("EnzoSimulationSerial::write");
}

//======================================================================

Timestep * 
EnzoSimulationSerial::create_timestep_ ( std::string name ) throw ()
/// @param name   Name of the timestep method to create (ignored)
{
  return new EnzoTimestep;
}

//----------------------------------------------------------------------

Initial * 
EnzoSimulationSerial::create_initial_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  
  Initial * initial = 0;

  if (name == "implosion_2d") {
    initial = new EnzoInitialImplosion2;
  }

  //--------------------------------------------------
  parameters_->set_current_group ("Initial");
  //--------------------------------------------------

  // parameter: Initial::cycle
  // parameter: Initial::time

  int cycle    = parameters_->value_integer ("cycle",0);
  double time  = parameters_->value_scalar ("time",0.0);

  ASSERT("EnzoSimulationSerial::create_initial_",
	 "create_mesh_ mush be called first",
	 mesh_ != NULL);

  ItPatch itPatch(mesh_);
  Patch * patch;
  while ((patch = ++itPatch)) {
    ItBlock itBlock (patch);
    Block * block;
    while ((block = ++itBlock)) {
      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);
      enzo_block->CycleNumber = cycle;
      enzo_block->Time        = time;
      enzo_block->OldTime     = time;
    }
  }

  return initial;
}

//----------------------------------------------------------------------

Stopping * 
EnzoSimulationSerial::create_stopping_ (std::string name) throw ()
/// @param name   Name of the stopping method to create (ignored)
{
  //--------------------------------------------------
  parameters_->set_current_group ("Stopping");
  //--------------------------------------------------

  // parameter: Stopping::cycle
  // parameter: Stopping::time

  int    stop_cycle = parameters_->value_integer
    ( "cycle" , std::numeric_limits<int>::max() );
  double stop_time  = parameters_->value_scalar
    ( "time" , std::numeric_limits<double>::max() );

  return new Stopping(stop_cycle,stop_time);
}

//----------------------------------------------------------------------

Boundary * 
EnzoSimulationSerial::create_boundary_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  
  Boundary * boundary = 0;

  return boundary;
}

//----------------------------------------------------------------------

Method * 
EnzoSimulationSerial::create_method_ ( std::string name ) throw ()
/// @param name   Name of the method to create
{

  Method * method = 0;

  if (name == "ppm")
    method = new EnzoMethodPpm  (parameters_);
  if (name == "ppml")
    method = new EnzoMethodPpml (parameters_);

  if (method == 0) {
    char buffer[80];
    sprintf (buffer,"Cannot create Method '%s'",name.c_str());
    ERROR("EnzoSimulationSerial::create_method", buffer);
  }

  return method;
}

//======================================================================

void EnzoSimulationSerial::output_images_
(
 Block * block,
 const char * file_format,
 int cycle,
 int cycle_skip
 ) throw ()
{

  if (! (cycle_skip && cycle % cycle_skip == 0)) return;

  Monitor * monitor = Monitor::instance();
  FieldBlock *       field_block = block->field_block();
  const FieldDescr * field_descr = field_block->field_descr();
  int nx,ny,nz;
  int gx,gy,gz;
  int mx,my,mz;
  field_block->enforce_boundary(boundary_reflecting);
  field_block->size(&nx,&ny,&nz);
  int count = field_descr->field_count();
  for (int index = 0; index < count; index++) {
    EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);
    field_descr->ghosts(index,&gx,&gy,&gz);
    mx=nx+2*gx;
    my=ny+2*gy;
    mz=nz+2*gz;
    char filename[255];
    std::string field_name = field_descr->field_name(index);
    Scalar * field_values = (Scalar *)field_block->field_values(index);
    sprintf (filename,file_format,
	     enzo_block->CycleNumber,index);
    monitor->image (filename, field_values, mx,my,mz, 2, reduce_sum, 
		    0.0, 1.0);
  }
}

//----------------------------------------------------------------------

void EnzoSimulationSerial::deallocate_() throw()
{
  delete stopping_;
  stopping_ = 0;
  delete timestep_;
  timestep_ = 0;
  delete initial_;
  initial_ = 0;
  for (size_t i=0; i<method_list_.size(); i++) {
    delete method_list_[i];
    method_list_[i] = 0;
  }
  
}
