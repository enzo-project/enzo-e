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

  if (initial == 0) {
    char buffer[80];
    sprintf (buffer,"Cannot create Initial '%s'",name.c_str());
    ERROR_MESSAGE("EnzoSimulation::create_initial", buffer);
  }

  return initial;
}

//----------------------------------------------------------------------

Boundary * 
EnzoSimulation::create_boundary_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  
  Boundary * boundary = 0;

  if (boundary == 0) {
    char buffer[80];
    sprintf (buffer,"Cannot create Boundary '%s'",name.c_str());
    ERROR_MESSAGE("EnzoSimulation::create_boundary", buffer);
  }

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

