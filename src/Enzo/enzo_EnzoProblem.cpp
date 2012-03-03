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

//======================================================================

