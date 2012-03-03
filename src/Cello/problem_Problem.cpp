// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Problem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "problem.hpp"

//----------------------------------------------------------------------

Problem::Problem() throw()
  : boundary_(0)
{
}

//----------------------------------------------------------------------

Problem::~Problem() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Problem::initialize_boundary(Parameters * parameters) throw()
{
  //--------------------------------------------------
  // parameter: Boundary : type
  //--------------------------------------------------

  std::string name = parameters->value_string("Boundary:type","");
  boundary_ = create_boundary_(name);
}

//----------------------------------------------------------------------

void Problem::deallocate_() throw()
{
  delete boundary_;      boundary_ = 0;
}

//----------------------------------------------------------------------

Boundary * Problem::create_boundary_ (std::string name) throw ()
{
  ERROR ("Problem::create_boundary_","Implictly abstract function called");
  return NULL;
}


