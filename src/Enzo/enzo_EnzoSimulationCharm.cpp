// $Id: enzo_EnzoSimulationCharm.cpp 2116 2011-03-17 22:11:13Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 Parameters * parameters,
 GroupProcess * group_process) throw ()
  : Simulation(new EnzoFactory,parameters,group_process)
{
}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::initialize() throw()
{
  // Call initialize for Simulation base class
  Simulation::initialize();

  // Initialize enzo namespace variables
  enzo::initialize(parameters_);

}

//----------------------------------------------------------------------

void EnzoSimulationCharm::finalize() throw()
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

    int    stop_reduce = reduce_patch->reduce_int (stop_block, reduce_op_land);

    stop_patch  = stop_patch && stop_reduce;
		      
    delete reduce_patch;

  }

  int    stop_mesh = reduce_mesh->reduce_int (stop_patch, reduce_op_land);
  
  //======================================================================
  // BEGIN MAIN LOOP
  //======================================================================

  while (! stop_mesh) {

    monitor->print("cycle %04d time %15.12f", cycle_,time_);

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

    ASSERT("EnzoSimulationCharm::run", "dt == 0", dt_mesh != 0.0);

    //--------------------------------------------------
    // Apply the methods
    //--------------------------------------------------

    // Initialize reductions

    int  stop_patch = true;

    while ((patch = ++it_patch)) {

      Reduce * reduce_patch = patch->group()->create_reduce ();

      int    stop_block  = true;

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

	// Initialize enzo_block aliases (may me modified)

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

  monitor->print("cycle %04d time %15.12f", cycle_,time_);

  performance.stop();

  performance.print();

}

//----------------------------------------------------------------------

void EnzoSimulationCharm::read() throw()
{
  INCOMPLETE("EnzoSimulationCharm::read");
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::write() const throw()
{
  INCOMPLETE("EnzoSimulationCharm::write");
}

//======================================================================

Timestep * 
EnzoSimulationCharm::create_timestep_ ( std::string name ) throw ()
/// @param name   Name of the timestep method to create (ignored)
{
  return new EnzoTimestep;
}

//----------------------------------------------------------------------

Initial * 
EnzoSimulationCharm::create_initial_ ( std::string name ) throw ()
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

  cycle_  = parameters_->value_integer ("cycle",0);
  time_   = parameters_->value_scalar ("time",0.0);

  ASSERT("EnzoSimulationCharm::create_initial_",
	 "create_mesh_ mush be called first",
	 mesh_ != NULL);

  ItPatch it_patch(mesh_);
  Patch * patch;
  while ((patch = ++it_patch)) {
    ItBlock it_block (patch);
    Block * block;
    while ((block = ++it_block)) {
      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);
      enzo_block->CycleNumber = cycle_;
      enzo_block->Time        = time_;
      enzo_block->OldTime     = time_;
    }
  }

  return initial;
}

//----------------------------------------------------------------------

Stopping * 
EnzoSimulationCharm::create_stopping_ (std::string name) throw ()
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
EnzoSimulationCharm::create_boundary_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  //--------------------------------------------------
  parameters_->set_current_group ("Boundary");
  //--------------------------------------------------

  // parameter: Boundary::type

  boundary_type_enum boundary_type = boundary_type_undefined;

  std::string boundary_param = parameters_->value_string("type", "undefined");

  if (       boundary_param ==   "reflecting") {
    boundary_type = boundary_type_reflecting;
  } else if (boundary_param ==   "outflow") {
    boundary_type = boundary_type_outflow;
  } else if (boundary_param ==   "inflow") {
    boundary_type = boundary_type_inflow;
  } else if (boundary_param ==   "periodic") {
    boundary_type = boundary_type_periodic;
  } else if (boundary_param ==   "file") {
    boundary_type = boundary_type_file;
  } else if (boundary_param ==   "code") {
    boundary_type = boundary_type_code;
  } else {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Illegal parameter Boundary::type '%s'",
	     boundary_param.c_str());
    ERROR("EnzoSimulationCharm::create_boundary_",
	  buffer);
  }
	     
  return new EnzoBoundary (boundary_type);
}

//----------------------------------------------------------------------

Method * 
EnzoSimulationCharm::create_method_ ( std::string name ) throw ()
/// @param name   Name of the method to create
{

  Method * method = 0;

  if (name == "ppm") {
    method = new EnzoMethodPpm  (parameters_);
  } else if (name == "ppml") {
    method = new EnzoMethodPpml (parameters_);
  }

  if (method == 0) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Cannot create Method '%s'",name.c_str());
    ERROR("EnzoSimulationCharm::create_method_", buffer);
  }

  return method;
}

//----------------------------------------------------------------------

Output * 
EnzoSimulationCharm::create_output_ ( std::string type ) throw ()
/// @param filename   File name format for the output object
{

  Output * output = 0;

  if (type == "image") {
    output = new EnzoOutputImage ();
  }

  if (output == 0) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Cannot create Output type '%s'",type.c_str());
    ERROR("EnzoSimulationCharm::create_output_", buffer);
  }

  return output;
}

