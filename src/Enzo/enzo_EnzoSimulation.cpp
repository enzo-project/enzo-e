// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation
(
 const char *   parameter_file,
#ifdef CONFIG_USE_CHARM
 int n,
#else
 GroupProcess * group_process,
#endif
 int            index
 ) throw()
  : Simulation
#ifdef CONFIG_USE_CHARM
    (parameter_file, n, index )
#else
    (parameter_file, group_process, index )
#endif
{
}

#ifdef CONFIG_USE_CHARM
/// Initialize an empty EnzoSimulation
EnzoSimulation::EnzoSimulation() 
{
  TRACE("EnzoSimulation()");
};

/// Initialize a migrated EnzoSimulation
EnzoSimulation::EnzoSimulation (CkMigrateMessage *m) 
{
  TRACE("EnzoSimulation(msg)");
};

//==================================================

#endif

//----------------------------------------------------------------------

EnzoSimulation::~EnzoSimulation() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize() throw()
{
  // Call initialize for Simulation base class
  Simulation::initialize();

  // Initialize enzo namespace variables
  enzo::initialize(parameters_);

}

//----------------------------------------------------------------------

void EnzoSimulation::finalize() throw()
{
}


//----------------------------------------------------------------------

void EnzoSimulation::run() throw()
{
  INCOMPLETE("EnzoSimulation::run");
}

//----------------------------------------------------------------------

void EnzoSimulation::read() throw()
{
  INCOMPLETE("EnzoSimulation::read");
}

//----------------------------------------------------------------------

void EnzoSimulation::write() const throw()
{
  INCOMPLETE("EnzoSimulation::write");
}

//----------------------------------------------------------------------

const Factory & EnzoSimulation::factory() const throw()
{ 
  if (factory_ == NULL) factory_ = new EnzoFactory;
  return *factory_;
}

//======================================================================

Stopping * EnzoSimulation::create_stopping_ (std::string name) throw ()
/// @param name   Name of the stopping method to create (ignored)
{
  //--------------------------------------------------
  // parameter: Stopping::cycle
  // parameter: Stopping::time
  //--------------------------------------------------

  int    stop_cycle = parameters_->value_integer
    ( "Stopping:cycle" , std::numeric_limits<int>::max() );
  double stop_time  = parameters_->value_float
    ( "Stopping:time" , std::numeric_limits<double>::max() );

  return new Stopping(stop_cycle,stop_time);
}

//----------------------------------------------------------------------

Timestep * EnzoSimulation::create_timestep_ ( std::string name ) throw ()
/// @param name   Name of the timestep method to create (ignored)
{
  return new EnzoTimestep;
}

//----------------------------------------------------------------------

Initial * EnzoSimulation::create_initial_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  
  Initial * initial = 0;

  if (name == "implosion_2d") {
    initial = new EnzoInitialImplosion2;
  }

  //--------------------------------------------------
  // parameter: Initial::cycle
  // parameter: Initial::time
  //--------------------------------------------------

  cycle_  = parameters_->value_integer ("Initial:cycle",0);
  time_   = parameters_->value_float   ("Initial:time",0.0);

  return initial;
}

//----------------------------------------------------------------------

Boundary * EnzoSimulation::create_boundary_ ( std::string name ) throw ()
/// @param name   Name of the initialization method to create
{
  //--------------------------------------------------
  // parameter: Boundary::type
  //--------------------------------------------------

  boundary_type_enum boundary_type = boundary_type_undefined;

  std::string boundary_param = 
    parameters_->value_string ("Boundary:type", "undefined");

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
    ERROR("EnzoSimulation::create_boundary_",
	  buffer);
  }
	     
  return new EnzoBoundary (boundary_type);
}

//----------------------------------------------------------------------

Method * EnzoSimulation::create_method_ ( std::string name ) throw ()
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
    ERROR("EnzoSimulation::create_method_", buffer);
  }

  return method;
}

//----------------------------------------------------------------------

Output * EnzoSimulation::create_output_ ( std::string type ) throw ()
/// @param filename   File name format for the output object
{

  Output * output = 0;

  if (type == "image") {
    output = new OutputImage ();
  }

  if (output == 0) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Cannot create Output type '%s'",type.c_str());
    ERROR("EnzoSimulation::create_output_", buffer);
  }

  return output;
}

