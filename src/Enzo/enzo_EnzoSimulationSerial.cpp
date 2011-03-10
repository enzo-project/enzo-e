// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationSerial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move boundary conditions to Boundary object; remove bc_reflecting
/// @todo     Move output to Output object
/// @todo     Move timestep reductions into Timestep object
/// @todo     Fix MPI reduction operations, which hang due to local block loops
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

  // @@@ performance_->start();
  Timer timer;
  Papi papi;
  
  timer.start();
  papi.start();

  //--------------------------------------------------
  // INITIALIZE FIELDS
  //--------------------------------------------------

  // Loop over Patches in the Mesh

  ItPatch itPatch(mesh_);
  EnzoPatch * patch;

  while ((patch = dynamic_cast <EnzoPatch*> (++itPatch))) {

    // Loop over Blocks in the Patch

    ItBlock itBlock(patch);
    EnzoBlock * enzo_block;

    while ((enzo_block = static_cast <EnzoBlock*> (++itBlock))) {

      // Initialize EnzoBlock field data

      initial_->compute(enzo_block);

      // Initialize EnzoBlock attributes

      enzo_block->initialize();

    }
  }

  //--------------------------------------------------
  // INITIAL STOPPING CRITERIA TEST
  //--------------------------------------------------

  int    cycle_patch = std::numeric_limits<int   >::max();
  double time_patch  = std::numeric_limits<double>::max();
  int    stop_patch  = true;

  while ((patch = dynamic_cast <EnzoPatch*> (++itPatch))) {

    ItBlock itBlock(patch);
    EnzoBlock * enzo_block;

    int    cycle_block = std::numeric_limits<int   >::max();
    double time_block  = std::numeric_limits<double>::max();
    int    stop_block  = true;

    while ((enzo_block = static_cast <EnzoBlock*> (++itBlock))) {
      cycle_block = MIN (cycle_block,enzo_block->CycleNumber);
      time_block  = MIN (time_block, enzo_block->Time);
      stop_block = stop_block && 
	stopping_->complete(enzo_block->CycleNumber,enzo_block->Time);
    }

#ifdef CONFIG_USE_MPI
    MPI_Allreduce (&cycle_block, &cycle_patch, 1, MPI_INT, MPI_MIN,
		   patch->layout()->mpi_comm());
    MPI_Allreduce (&time_block, &time_patch, 1, MPI_DOUBLE, MPI_MIN,
		   patch->layout()->mpi_comm());
    MPI_Allreduce (&stop_block, &stop_patch, 1, MPI_INT, MPI_LAND,
		   patch->layout()->mpi_comm());
#else
    cycle_patch = cycle_block;
    time_patch = time_block;
    stop_patch = stop_block;
#endif

  }

  int    cycle_mesh;
  double time_mesh;
  int    stop_mesh;
  
#ifdef CONFIG_USE_MPI
    MPI_Allreduce (&cycle_patch, &cycle_mesh, 1, MPI_INT, MPI_MIN,
		   mesh_->mpi_comm());
    MPI_Allreduce (&time_patch, &time_mesh, 1, MPI_DOUBLE, MPI_MIN,
		   mesh_->mpi_comm());
    MPI_Allreduce (&stop_patch, &stop_mesh, 1, MPI_INT, MPI_LAND,
		   mesh_->mpi_comm());
#else
  cycle_mesh = cycle_patch;
  time_mesh = time_patch;
  stop_mesh = stop_patch;
#endif

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

    while ((patch = dynamic_cast <EnzoPatch*> (++itPatch))) {

      ItBlock itBlock(patch);
      EnzoBlock * enzo_block;

      double dt_block = std::numeric_limits<double>::max();

      // Accumulate Block-local timesteps

      while ((enzo_block = static_cast <EnzoBlock*> (++itBlock))) {

	// Enforce boundary conditions
	enzo_block->field_block()->enforce_boundary(boundary_reflecting);

	// Output while we're here
	output_images_(enzo_block, "enzo-p-%06d.%d.png",cycle_mesh,10);

	// Accumulate Block-local dt
	dt_block = MIN(dt_block,timestep_->compute(enzo_block));

      } // ( enzo_block = ++itBlock )

      // Reduce Block timesteps to Patch

#ifdef CONFIG_USE_MPI
      MPI_Allreduce (&dt_block, &dt_patch, 1, MPI_DOUBLE, MPI_MIN,
		     patch->layout()->mpi_comm());
#else
      dt_patch = dt_block;
#endif

    } // ( patch = ++itPatch )

    // Reduce Patch timesteps to Mesh

    double dt_mesh;

#ifdef CONFIG_USE_MPI
    MPI_Allreduce (&dt_patch, &dt_mesh, 1, MPI_DOUBLE, MPI_MIN,
		   mesh_->mpi_comm());
#else
    dt_mesh = dt_patch;
#endif

    ASSERT("EnzoSimulationSerial::run", "dt == 0", dt_mesh != 0.0);

    // Assign the computed timestep

    //--------------------------------------------------
    // Apply the methods
    //--------------------------------------------------

    // Initialize reductions

    cycle_patch = std::numeric_limits<int   >::max();
    time_patch  = std::numeric_limits<double>::max();
    stop_patch  = true;

    while ((patch = dynamic_cast <EnzoPatch*> (++itPatch))) {

      ItBlock itBlock(patch);
      EnzoBlock * enzo_block;

      int    cycle_block = std::numeric_limits<int   >::max();
      double time_block  = std::numeric_limits<double>::max();
      int    stop_block  = true;

      while ((enzo_block = static_cast <EnzoBlock*> (++itBlock))) {

	// UNIFORM TIMESTEP OVER ALL BLOCKS IN MESH

	enzo_block->dt = MIN(dt_mesh, (stopping_->stop_time() - enzo_block->dt));
	
	// Loop through methods

	for (size_t i = 0; i < method_list_.size(); i++) {
	  Method * method = method_list_[i];

	  method -> compute_block (enzo_block,enzo_block->Time,enzo_block->dt);

	} // method

	// Update cycle and time

	enzo_block->CycleNumber += 1;
	enzo_block->OldTime      = enzo_block->Time;
	enzo_block->Time        += enzo_block->dt;

	// Global cycle and time reduction

	cycle_block = MIN (cycle_block,enzo_block->CycleNumber);
	time_block  = MIN (time_block, enzo_block->Time);
	stop_block = stop_block && 
	  stopping_->complete(enzo_block->CycleNumber,enzo_block->Time);

      } // (enzo_block = ++itBlock)
#ifdef CONFIG_USE_MPI
      MPI_Allreduce (&cycle_block, &cycle_patch, 1, MPI_INT, MPI_MIN,
		     patch->layout()->mpi_comm());
      MPI_Allreduce (&time_block, &time_patch, 1, MPI_DOUBLE, MPI_MIN,
		     patch->layout()->mpi_comm());
      MPI_Allreduce (&stop_block, &stop_patch, 1, MPI_INT, MPI_LAND,
		     patch->layout()->mpi_comm());
#else
      cycle_patch = cycle_block;
      time_patch = time_block;
      stop_patch = stop_block;
#endif
    } // (patch = ++itPatch)
#ifdef CONFIG_USE_MPI
    MPI_Allreduce (&cycle_patch, &cycle_mesh, 1, MPI_INT, MPI_MIN,
		   mesh_->mpi_comm());
    MPI_Allreduce (&time_patch, &time_mesh, 1, MPI_DOUBLE, MPI_MIN,
		   mesh_->mpi_comm());
    MPI_Allreduce (&stop_patch, &stop_mesh, 1, MPI_INT, MPI_LAND,
		   mesh_->mpi_comm());
#else
    cycle_mesh = cycle_patch;
    time_mesh = time_patch;
    stop_mesh = stop_patch;
#endif

  } // ! stop

  //======================================================================
  // END MAIN LOOP
  //======================================================================

  monitor->print("cycle %04d time %15.12f", cycle_mesh,time_mesh);

  while ((patch = dynamic_cast <EnzoPatch*> (++itPatch))) {

    ItBlock itBlock(patch);
    EnzoBlock * enzo_block;

    while ((enzo_block = static_cast <EnzoBlock*> (++itBlock))) {

      output_images_(enzo_block, "enzo-p-%06d.%d.png",cycle_mesh,1);

    }
  }

  // @@@ performance_->stop();
  papi.stop();
  timer.stop();

#ifdef CONFIG_USE_PAPI
  monitor->print ("PAPI Time real   = %f",papi.time_real());
  monitor->print ("PAPI Time proc   = %f\n",papi.time_proc());
  monitor->print ("PAPI GFlop count = %f\n",papi.flop_count()*1e-9);
  //  MONITOR->PRINT ("PAPI GFlop rate  = %f\n",papi.flop_rate()*1e-9);
  monitor->print ("PAPI GFlop rate  = %f\n",
		   1e-9 * papi.flop_count() / papi.time_real());
  
#endif

  monitor->print ("Real time = %f\n",timer.value());

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
  EnzoPatch * patch;
  while ((patch = dynamic_cast <EnzoPatch*> (++itPatch))) {
    ItBlock itBlock (patch);
    EnzoBlock * enzo_block;
    while ((enzo_block = static_cast <EnzoBlock*> (++itBlock))) {
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

void EnzoSimulationSerial::block_start_ (EnzoBlock * enzo_block) throw ()
{

}

//----------------------------------------------------------------------

void EnzoSimulationSerial::output_images_
(
 EnzoBlock * enzo_block,
 const char * file_format,
 int cycle,
 int cycle_skip
 ) throw ()
{

  if (! (cycle_skip && cycle % cycle_skip == 0)) return;

  Monitor * monitor = Monitor::instance();
  FieldBlock *       field_block = enzo_block->field_block();
  const FieldDescr * field_descr = field_block->field_descr();
  int nx,ny,nz;
  int gx,gy,gz;
  int mx,my,mz;
  field_block->enforce_boundary(boundary_reflecting);
  field_block->size(&nx,&ny,&nz);
  int count = field_descr->field_count();
  for (int index = 0; index < count; index++) {
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
