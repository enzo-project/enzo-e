// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation(Error   * error,
			       Monitor * monitor) throw ()
  : Simulation(error,monitor),
    enzo_descr_(new EnzoDescr())
{
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize(FILE * fp) throw()
{
  // Call initialize for Simulation base class
  Simulation::initialize(fp);

  // Call initialize for Enzo-specific Simulation
  enzo_descr_->initialize(parameters_);

  // @@@ WRITE OUT ENZO DESCRIPTION FOR DEBUGGING
  enzo_descr_->write(stdout);
}

//----------------------------------------------------------------------

void EnzoSimulation::finalize() throw()
{
  delete enzo_descr_;
}

//----------------------------------------------------------------------

void EnzoSimulation::run() throw()
{
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  //   FileHdf5   hdf5;

  // output_progress(monitor,cycle,time,dt);
  // output_images(monitor,cycle,field_block);

  // Papi papi;
  // Timer timer;
  
  // timer.start();
  // papi.start();

  control_->initialize(data_descr());

  while (! control_->is_done() ) {

    {
      char buffer[40];
      sprintf (buffer,"time %g cycle %d\n",
	       enzo_descr_->Time,
	       enzo_descr_->CycleNumber);
      monitor_->print(buffer);
    }

    control_->initialize_cycle();

    // for each Block in Patch

    ItBlocks itBlocks(mesh_->root_patch());
    while (DataBlock * data_block = (DataBlock *) ++itBlocks) {

      control_->initialize_block(data_block);
      
      control_->finalize_block(data_block);

    } // Block in Patch

    control_->finalize_cycle();
  }

  control_->finalize(data_descr());

  // while (time < time_final && cycle <= cycle_final) {

  //   simulation->initialize_block(data_block);

  //   field_block->enforce_boundary(boundary_reflecting);

  //   dt = simulation->timestep()->compute_block(data_block);

  //   output_images  (monitor,cycle,field_block,cycle_image);
  //   output_progress(monitor,cycle,time,dt,cycle_progress);
  //   output_dump    (hdf5,cycle,field_block,cycle_dump);

  //   simulation->method(0)->compute_block(data_block,time,dt);

  //   simulation->finalize_block(data_block);

  //   ++cycle;
  //   time += MIN(dt,time_final-time);

  // }

  // papi.stop();
  // timer.stop();

  // output_images  (monitor,cycle_final,field_block);
  // output_progress(monitor,cycle_final,time,dt);

  // printf ("Time real = %f\n",papi.time_real());
  // printf ("Time proc = %f\n",papi.time_proc());
  // printf ("Flop count = %lld\n",papi.flop_count());
  // printf ("GFlop rate = %f\n",papi.flop_rate()*1e-9);

  // printf ("Timer time = %f\n",timer.value());

}

//----------------------------------------------------------------------

void EnzoSimulation::read() throw()
{
  INCOMPLETE_MESSAGE("EnzoSimulation::read","");
}

//----------------------------------------------------------------------

void EnzoSimulation::write() throw()
{
  INCOMPLETE_MESSAGE("EnzoSimulation::write","");
}

//======================================================================

Control * 
EnzoSimulation::create_control_ (std::string name) throw ()
/// @param name   Name of the control method to create (ignored)
{
  return new EnzoControl(error_, monitor_,parameters_,enzo_descr_);
}

//----------------------------------------------------------------------

Timestep * 
EnzoSimulation::create_timestep_ ( std::string name ) throw ()
/// @param name   Name of the timestep method to create (ignored)
{
  return new EnzoTimestep(enzo_descr_);
}

//----------------------------------------------------------------------

Initial * 
EnzoSimulation::create_initial_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  
  Initial * initial = 0;

  if (name == "implosion2")  
    initial = new EnzoInitialImplosion2 (error_, monitor_, enzo_descr_);

  return initial;
}

//----------------------------------------------------------------------

Boundary * 
EnzoSimulation::create_boundary_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  
  Boundary * boundary = 0;

  return boundary;
}

//----------------------------------------------------------------------

Method * 
EnzoSimulation::create_method_ ( std::string name ) throw ()
/// @param name   Name of the method to create
{

  Method * method = 0;

  if (name == "ppm")
    method = new EnzoMethodPpm  (error_,monitor_,parameters_,enzo_descr_);
  if (name == "ppml")
    method = new EnzoMethodPpml (error_,monitor_,parameters_,enzo_descr_);

  if (method == 0) {
    char buffer[80];
    sprintf (buffer,"Cannot create Method '%s'",name.c_str());
    ERROR_MESSAGE("EnzoSimulation::create_method", buffer);
  }

  return method;
}

