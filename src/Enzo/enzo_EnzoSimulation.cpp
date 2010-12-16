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
    enzo_(new EnzoDescr(global))
{
  initialize_enzo();
}

//======================================================================

void EnzoSimulation::initialize_enzo () throw ()
{
  Parameters * parameters = global_->parameters();

  // --------------------------------------------------
  // Initiazize user descriptors
  // --------------------------------------------------

  // --------------------------------------------------
  // Initiazize user method descriptors
  // --------------------------------------------------

  parameters->set_current_group("Method");

  int method_count = parameters->list_length("sequence");

  if (method_count == 0) {
    ERROR_MESSAGE ("Simulation::initialize",
		   "List parameter 'Method sequence' must have length greater than zero");
  }

  for (int i=0; i<method_count; i++) {

    std::string method_name = parameters->list_value_string(i,"sequence");

    MethodHyperbolic * method_hyperbolic = add_method(method_name);

    method_hyperbolic->initialize(data_);
  }


  // --------------------------------------------------
  // Initialize method control
  // --------------------------------------------------

  set_control("ignored");

  // --------------------------------------------------
  // Initialize method timestep
  // --------------------------------------------------

  set_timestep("ignored");

}

//----------------------------------------------------------------------

MethodControl * 
EnzoSimulation::create_control_ (std::string control_name) throw ()
/// @param control_name   Name of the control method to create
{
  return new MethodEnzoControl(global_,enzo_);
}

//----------------------------------------------------------------------

MethodTimestep * 
EnzoSimulation::create_timestep_ ( std::string timestep_name ) throw ()
/// @param timestep_name   Name of the timestep method to create
{
  return new MethodEnzoTimestep(enzo_);
}

//----------------------------------------------------------------------

MethodHyperbolic * 
EnzoSimulation::create_method_ ( std::string method_name ) throw ()
/// @param method_name   Name of the method to create
{
  printf ("%s:%d '%s'\n",__FILE__,__LINE__,method_name.c_str());
  if (method_name == "ppm") {
    return new MethodEnzoPpm (global_,enzo_);
  } else {
    char buffer[80];
    sprintf (buffer,"Unknown method '%s'",method_name.c_str());
    global_->error()->error_ (__FILE__,__LINE__,"EnzoSimulation::create_method",
				       buffer);
    return 0;
  }
}

