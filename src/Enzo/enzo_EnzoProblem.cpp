// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProblem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-03
/// @brief    Implementation of EnzoProblem class
///
/// 

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoProblem::EnzoProblem() throw ()
{
  INCOMPLETE("EnzoProblem::EnzoProblem");
}

//----------------------------------------------------------------------

EnzoProblem::~EnzoProblem() throw ()
{
  INCOMPLETE("EnzoProblem::~EnzoProblem");
}

//----------------------------------------------------------------------

EnzoProblem::EnzoProblem(const EnzoProblem & EnzoProblem) throw ()
/// @param     EnzoProblem  Object being copied
{
  INCOMPLETE("EnzoProblem::EnzoProblem(EnzoProblem)");
}

//----------------------------------------------------------------------

EnzoProblem & EnzoProblem::operator= (const EnzoProblem & EnzoProblem) throw ()
/// @param     EnzoProblem  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE("EnzoProblem::operator=");
  return *this;
}

//======================================================================

Boundary * EnzoProblem::create_boundary_ ( std::string name ) throw ()
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
    ERROR1("EnzoProblem::create_boundary_",
	   "Unrecognized boundary type '%s'",
	   name.c_str());
  }
	     
  return new EnzoBoundary (boundary_type);
}

//----------------------------------------------------------------------

Stopping * EnzoProblem::create_stopping_ 
(
 std::string name,
 int         stop_cycle,
 double      stop_time
 ) throw ()
/// @param name   Name of the stopping method to create (ignored)
/// @param stop_cycle  Stopping cycle
/// @param stop_time  Stopping time
{

  return new Stopping(stop_cycle,stop_time);
}

//----------------------------------------------------------------------

Timestep * EnzoProblem::create_timestep_ ( std::string name ) throw ()
/// @param name   Name of the timestep method to create (ignored)
{
  if (name == "ppml") {
    return new EnzoTimestepPpml;
  } else {
    return new EnzoTimestep;
  }
}

//----------------------------------------------------------------------

Initial * EnzoProblem::create_initial_ 
(
 std::string name,
 int         init_cycle,
 double      init_time
 ) throw ()
/// @param name   Name of the initialization method to create
{
  
  Initial * initial = 0;

  if (name == "implosion_2d") {
    initial = new EnzoInitialImplosion2(init_cycle,init_time);
  }

  return initial;
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------

Method * EnzoProblem::create_method_ 
(
 std::string name,
 Parameters * parameters) throw ()
/// @param name   Name of the method to create
{

  Method * method = 0;

  if (name == "ppm") {
    method = new EnzoMethodPpm  (parameters);
  } else if (name == "ppml") {
    method = new EnzoMethodPpml (parameters);
  }

  if (method == 0) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Cannot create Method '%s'",name.c_str());
    ERROR("EnzoProblem::create_method_", buffer);
  }

  return method;
}

//----------------------------------------------------------------------

Output * EnzoProblem::create_output_ 
(
 std::string    type,
 GroupProcess * group_process,
 Hierarchy    * hierarchy,
 const Factory * factory
 ) throw ()
/// @param type          Type of Output object to create
/// @param group_process Image output needs group process size
/// @param hierarchy     Image output needs image size (currently patch(0) size)
/// @param filename   File name format for the output object
{

  Output * output = NULL;

  // Create Enzo-specific output types

  // ...

  // Create a Cello output type using base class

  if (output == NULL) {
    output = Problem::create_output_(type,group_process,hierarchy,factory);
  }

  if (output == NULL) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Cannot create Output type '%s'",type.c_str());
    ERROR("EnzoProblem::create_output_", buffer);
  }

  return output;
}

