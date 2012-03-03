// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
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
 CProxy_BlockReduce proxy_block_reduce
#else
 GroupProcess * group_process
#endif
 ) throw()
  : Simulation
#ifdef CONFIG_USE_CHARM
    (parameter_file, n, proxy_block_reduce)
#else
    (parameter_file, group_process)
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
  enzo::initialize(parameters_,field_descr());

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
  // parameter: Stopping : cycle
  // parameter: Stopping : time
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
  if (name == "ppml") {
    return new EnzoTimestepPpml;
  } else {
    return new EnzoTimestep;
  }
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
  // parameter: Initial : cycle
  //--------------------------------------------------

  cycle_  = parameters_->value_integer ("Initial:cycle",0);

  //--------------------------------------------------
  // parameter: Initial : time
  //--------------------------------------------------

  time_   = parameters_->value_float   ("Initial:time",0.0);

  return initial;
}

//----------------------------------------------------------------------

Boundary * EnzoSimulation::create_boundary_ ( std::string name ) throw ()
/// @param name   Name of boundary condition to use
{

  boundary_type_enum boundary_type = boundary_type_undefined;

  if (       name == "reflecting") {
    boundary_type = boundary_type_reflecting;
  } else if (name == "outflow") {
    boundary_type = boundary_type_outflow;
  } else if (name == "inflow") {
    boundary_type = boundary_type_inflow;
  } else if (name == "periodic") {
    boundary_type = boundary_type_periodic;
  } else {
    ERROR1("EnzoSimulation::create_boundary_",
	   "Unrecognized boundary type '%s'",
	   name.c_str());
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

  Output * output = NULL;

  // Create Enzo-specific output types

  // ...

  // Create a Cello output type using base class

  if (output == NULL) {
    output = Simulation::create_output_(type);
  }

  if (output == NULL) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Cannot create Output type '%s'",type.c_str());
    ERROR("EnzoSimulation::create_output_", buffer);
  }

  return output;
}

