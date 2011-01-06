// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

#include "cello.hpp"
#include "enzo.hpp"

#include "parameters.hpp"
#include "simulation.hpp"

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation(Global * global) throw ()
  : Simulation(global),
    enzo_descr_(0)
{
  TRACE_MESSAGE;
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize(std::string parameter_file) throw()
{
  // Call initialize for Simulation base class

  Simulation::initialize(parameter_file);

  // Call initialize for Enzo-specific Simulation
  initialize_enzo_();
}

//======================================================================

void EnzoSimulation::initialize_enzo_() throw ()
{
  enzo_descr_ = new EnzoDescr(global_);
}

//----------------------------------------------------------------------

MethodControl * 
EnzoSimulation::create_control_ (std::string control_name) throw ()
/// @param control_name   Name of the control method to create
{
  return new MethodEnzoControl(global_,enzo_descr_);
}

//----------------------------------------------------------------------

MethodTimestep * 
EnzoSimulation::create_timestep_ ( std::string timestep_name ) throw ()
/// @param timestep_name   Name of the timestep method to create
{
  return new MethodEnzoTimestep(enzo_descr_);
}

//----------------------------------------------------------------------

MethodHyperbolic * 
EnzoSimulation::create_method_ ( std::string method_name ) throw ()
/// @param method_name   Name of the method to create
{

  MethodHyperbolic * method = 0;

  if (method_name == "ppm")
    method = new MethodEnzoPpm  (global_,enzo_descr_);
  if (method_name == "ppml")
    method = new MethodEnzoPpml (global_,enzo_descr_);

  if (method == 0) {
    char buffer[80];
    sprintf (buffer,"Cannot create Method '%s'",method_name.c_str());
    ERROR_MESSAGE("EnzoSimulation::create_method", buffer);
  }

  return method;
}

//----------------------------------------------------------------------

MethodInitial * 
EnzoSimulation::create_initial_ ( std::string initial_name ) throw ()
/// @param initial_name   Name of the initialization method to create
{
  
  MethodInitial * initial = 0;

  if (initial_name == "implosion2")  
    initial = new MethodInitialImplosion2 (global_,enzo_descr_);

  if (initial == 0) {
    char buffer[80];
    sprintf (buffer,"Cannot create Initialization '%s'",initial_name.c_str());
    ERROR_MESSAGE("EnzoSimulation::create_initial", buffer);
  }

  return initial;
}

