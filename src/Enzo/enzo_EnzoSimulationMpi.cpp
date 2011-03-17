// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulationMpi user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationMpi::EnzoSimulationMpi
(
 Parameters * parameters,
 GroupProcess * group_process) throw ()
  : Simulation(new EnzoFactory,parameters,group_process)
{
}

//----------------------------------------------------------------------

EnzoSimulationMpi::~EnzoSimulationMpi() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::initialize() throw()
{
  // Call initialize for Simulation base class
  Simulation::initialize();

  // Initialize enzo namespace variables
  enzo::initialize(parameters_);

}

//----------------------------------------------------------------------

void EnzoSimulationMpi::finalize() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::run() throw()
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

  //--------------------------------------------------
  // INITIAL STOPPING CRITERIA TEST
  //--------------------------------------------------

  Reduce * reduce_mesh = mesh_->group()->create_reduce ();

  int    cycle_patch = std::numeric_limits<int>::max();
  int    stop_patch  = true;
  double time_patch  = std::numeric_limits<double>::max();

  while ((patch = ++it_patch)) {

    Reduce * reduce_patch = patch->group()->create_reduce ();

    int    cycle_block = std::numeric_limits<int>::max();
    int    stop_block  = true;
    double time_block  = std::numeric_limits<double>::max();

    ItBlock it_block(patch);
    Block * block;

    while ((block = ++it_block)) {

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

    ASSERT("EnzoSimulationMpi::run", "dt == 0", dt_mesh != 0.0);

    //--------------------------------------------------
    // Apply the methods
    //--------------------------------------------------

    // Initialize reductions

    int    cycle_patch = std::numeric_limits<int>::max();
    int    stop_patch  = true;
    double time_patch  = std::numeric_limits<double>::max();

    while ((patch = ++it_patch)) {

      Reduce * reduce_patch = patch->group()->create_reduce ();

      int    cycle_block = std::numeric_limits<int>::max();
      int    stop_block  = true;
      enzo_float time_block  = std::numeric_limits<enzo_float>::max();

      ItBlock it_block(patch);
      Block * block;

      while ((block = ++it_block)) {

	EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

	// Initialize enzo_block aliases (may me modified)

	int &    cycle        = enzo_block->CycleNumber;
	enzo_float & time     = enzo_block->Time;
	enzo_float & dt       = enzo_block->dt;
	enzo_float & old_time = enzo_block->OldTime;

	// Perform any scheduled output

	for (size_t i=0; i<output_list_.size(); i++) {
	  Output * output = output_list_[i];
	  output->scheduled_write(block,cycle,time);
	}

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

      } // (block = ++it_block)

      int cycle_reduce = reduce_patch->reduce_int (cycle_block, reduce_op_min);
      double time_reduce = reduce_patch->reduce_double (time_block, reduce_op_min);
      int    stop_reduce = reduce_patch->reduce_int (stop_block, reduce_op_land);

      cycle_patch = MIN(cycle_patch, cycle_reduce);
      time_patch  = MIN(time_patch, time_reduce);
      stop_patch  = stop_patch && stop_reduce;

      delete reduce_patch;
    } // (patch = ++it_patch)

    cycle_mesh = reduce_mesh->reduce_int (cycle_patch, reduce_op_min);
    time_mesh  = reduce_mesh->reduce_double (time_patch, reduce_op_min);
    stop_mesh  = reduce_mesh->reduce_int (stop_patch, reduce_op_land);

  } // ! stop

  //======================================================================
  // END MAIN LOOP
  //======================================================================

  monitor->print("cycle %04d time %15.12f", cycle_mesh,time_mesh);

  performance.stop();

  performance.print();

}

//----------------------------------------------------------------------

void EnzoSimulationMpi::read() throw()
{
  INCOMPLETE("EnzoSimulationMpi::read");
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::write() const throw()
{
  INCOMPLETE("EnzoSimulationMpi::write");
}

//======================================================================

Timestep * 
EnzoSimulationMpi::create_timestep_ ( std::string name ) throw ()
/// @param name   Name of the timestep method to create (ignored)
{
  return new EnzoTimestep;
}

//----------------------------------------------------------------------

Initial * 
EnzoSimulationMpi::create_initial_ ( std::string name ) throw ()
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

  ASSERT("EnzoSimulationMpi::create_initial_",
	 "create_mesh_ mush be called first",
	 mesh_ != NULL);

  ItPatch it_patch(mesh_);
  Patch * patch;
  while ((patch = ++it_patch)) {
    ItBlock it_block (patch);
    Block * block;
    while ((block = ++it_block)) {
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
EnzoSimulationMpi::create_stopping_ (std::string name) throw ()
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
EnzoSimulationMpi::create_boundary_ ( std::string name ) throw ()
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
    ERROR("EnzoSimulationMpi::create_boundary_",
	  buffer);
  }
	     
  return new EnzoBoundary (boundary_type);
}

//----------------------------------------------------------------------

Method * 
EnzoSimulationMpi::create_method_ ( std::string name ) throw ()
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
    ERROR("EnzoSimulationMpi::create_method_", buffer);
  }

  return method;
}

//----------------------------------------------------------------------

Output * 
EnzoSimulationMpi::create_output_ ( std::string type ) throw ()
/// @param filename   File name format for the output object
{

  Output * output = 0;

  if (type == "image") {
    output = new EnzoOutputImage ();
  }

  if (output == 0) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Cannot create Output type '%s'",type.c_str());
    ERROR("EnzoSimulationMpi::create_output_", buffer);
  }

  return output;
}

//======================================================================

void EnzoSimulationMpi::output_images_
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

void EnzoSimulationMpi::deallocate_() throw()
{
  
}
