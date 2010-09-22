// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

#include "cello.hpp"
#include "global.hpp"
#include "parameters.hpp"
#include "simulation.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation(Global * global) throw ()
  : Simulation(global),
    enzo_(new EnzoDescr(global))
{
  user_initialize_();
}

//======================================================================

void EnzoSimulation::user_initialize_ () throw ()
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

    UserMethod * user_method = add_user_method(method_name);

    user_method->initialize(data_descr_);
  }


  // --------------------------------------------------
  // Initialize user control
  // --------------------------------------------------

  set_user_control("ignored");

  // --------------------------------------------------
  // Initialize user timestep
  // --------------------------------------------------

  set_user_timestep("ignored");

}

//----------------------------------------------------------------------

UserControl * 
EnzoSimulation::create_user_control_ (std::string control_name) throw ()
{
  return new MethodEnzoControl(global_,enzo_);
}

//----------------------------------------------------------------------

UserTimestep * 
EnzoSimulation::create_user_timestep_ ( std::string timestep_name ) throw ()
{
  return new MethodEnzoTimestep(enzo_);
}

//----------------------------------------------------------------------

UserMethod * 
EnzoSimulation::create_user_method_ ( std::string method_name ) throw ()
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

