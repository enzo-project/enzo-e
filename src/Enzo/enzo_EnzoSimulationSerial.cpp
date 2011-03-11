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

EnzoSimulationSerial::EnzoSimulationSerial(Parameters * parameters) throw ()
  : Simulation(parameters)
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

    }
  }

  //--------------------------------------------------
  // INITIAL STOPPING CRITERIA TEST
  //--------------------------------------------------

#ifdef CONFIG_USE_MPI
  Reduce<int>       reduce_cycle_patch(op_min,MPI_INT, mesh_->mpi_comm());
  Reduce<double> reduce_time_patch (op_min,MPI_DOUBLE, mesh_->mpi_comm());
  Reduce<int>       reduce_stop_patch (op_land,MPI_INT,mesh_->mpi_comm());
#else
  Reduce<int>    reduce_cycle_patch(op_min);
  Reduce<double> reduce_time_patch (op_min);
  Reduce<int>    reduce_stop_patch (op_land);
#endif

  while ((patch = ++itPatch)) {

#ifdef CONFIG_USE_MPI
    Reduce<int>       reduce_cycle_block(op_min,MPI_INT, patch->layout()->mpi_comm());
    Reduce<double> reduce_time_block (op_min,MPI_DOUBLE, patch->layout()->mpi_comm());
    Reduce<int>       reduce_stop_block (op_land,MPI_INT,patch->layout()->mpi_comm());
#else
    Reduce<int>    reduce_cycle_block(op_min);
    Reduce<double> reduce_time_block (op_min);
    Reduce<int>    reduce_stop_block (op_land);
#endif

    ItBlock itBlock(patch);
    Block * block;

    while ((block = ++itBlock)) {

      EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

      reduce_cycle_block.accum(enzo_block->CycleNumber);
      reduce_time_block.accum (enzo_block->Time);
      reduce_stop_block.accum 
	(stopping_->complete(enzo_block->CycleNumber,enzo_block->Time));
    }

    reduce_time_patch.accum  (reduce_time_block.reduce());
    reduce_cycle_patch.accum (reduce_cycle_block.reduce());
    reduce_stop_patch.accum  (reduce_stop_block.reduce());
    
  }

  int    cycle_mesh = reduce_cycle_patch.reduce();
  double time_mesh  = reduce_time_patch.reduce();
  int    stop_mesh  = reduce_stop_patch.reduce();

  
  //======================================================================
  // BEGIN MAIN LOOP
  //======================================================================

  while (! stop_mesh) {

    monitor->print("cycle %04d time %15.12f", cycle_mesh,time_mesh);

    //--------------------------------------------------
    // Determine timestep and dump output
    //--------------------------------------------------

#ifdef CONFIG_USE_MPI
    Reduce<double> reduce_dt_patch (op_min,MPI_DOUBLE,mesh_->mpi_comm());
#else
    Reduce<double> reduce_dt_patch (op_min);
#endif

    // Accumulate Patch-local timesteps

    while ((patch = ++itPatch)) {

#ifdef CONFIG_USE_MPI
      Reduce<double> reduce_dt_block(op_min, MPI_DOUBLE, patch->layout()->mpi_comm());
#else
      Reduce<double> reduce_dt_block(op_min);
#endif

      ItBlock itBlock(patch);
      Block * block;

      // Accumulate Block-local timesteps

      while ((block = ++itBlock)) {

	// Enforce boundary conditions
	block->field_block()->enforce_boundary(boundary_reflecting);

	// Output while we're here
	output_images_(block, "enzo-p-%06d.%d.png",cycle_mesh,0);

	// Accumulate Block-local dt

	reduce_dt_block.accum (timestep_->compute(block));

      } // ( block = ++itBlock )

      reduce_dt_patch.accum(reduce_dt_block.reduce());

    } // ( patch = ++itPatch )

    double dt_mesh  = reduce_dt_patch.reduce();

    ASSERT("EnzoSimulationSerial::run", "dt == 0", dt_mesh != 0.0);

    //--------------------------------------------------
    // Apply the methods
    //--------------------------------------------------

    // Initialize reductions

    reduce_cycle_patch.reset();
    reduce_time_patch.reset();
    reduce_stop_patch.reset();

    while ((patch = ++itPatch)) {

#ifdef CONFIG_USE_MPI
      Reduce<int>    reduce_cycle_block(op_min,MPI_INT, patch->layout()->mpi_comm());
      Reduce<double> reduce_time_block (op_min,MPI_DOUBLE, patch->layout()->mpi_comm());
      Reduce<int>    reduce_stop_block (op_land,MPI_INT,patch->layout()->mpi_comm());
#else
      Reduce<int>    reduce_cycle_block(op_min);
      Reduce<double> reduce_time_block (op_min);
      Reduce<int>    reduce_stop_block (op_land);
#endif

      ItBlock itBlock(patch);
      Block * block;

      while ((block = ++itBlock)) {

	EnzoBlock * enzo_block = static_cast <EnzoBlock*> (block);

	// UNIFORM TIMESTEP OVER ALL BLOCKS IN MESH

	double dt_stop = (stopping_->stop_time() - enzo_block->dt);

	enzo_block->dt = MIN(dt_mesh, dt_stop);
	
	// Loop through methods

	for (size_t i = 0; i < method_list_.size(); i++) {

	  Method * method = method_list_[i];

	  method -> compute_block (block,enzo_block->Time,enzo_block->dt);

	} // method

	// Update cycle and time

	enzo_block->CycleNumber += 1;
	enzo_block->OldTime      = enzo_block->Time;
	enzo_block->Time        += enzo_block->dt;

	// Global cycle and time reduction

	reduce_cycle_block.accum(enzo_block->CycleNumber);
	reduce_time_block.accum (enzo_block->Time);
	reduce_stop_block.accum 
	  (stopping_->complete(enzo_block->CycleNumber,enzo_block->Time));


      } // (block = ++itBlock)

      reduce_cycle_patch.accum (reduce_cycle_block.reduce());
      reduce_time_patch.accum  (reduce_time_block.reduce());
      reduce_stop_patch.accum  (reduce_stop_block.reduce());


    } // (patch = ++itPatch)

    cycle_mesh = reduce_cycle_patch.reduce();
    time_mesh  = reduce_time_patch.reduce();
    stop_mesh  = reduce_stop_patch.reduce();

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
  INCOMPLETE("EnzoSimulationSerial::read","");
}

//----------------------------------------------------------------------

void EnzoSimulationSerial::write() throw()
{
  INCOMPLETE("EnzoSimulationSerial::write","");
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
