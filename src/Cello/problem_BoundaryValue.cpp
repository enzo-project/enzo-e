// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryValue.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-02
/// @brief    Implementation of the default BoundaryValue boundary value class

#include "problem.hpp"

//----------------------------------------------------------------------

void BoundaryValue::enforce 
(const FieldDescr * field_descr,
 CommBlock * block,
 face_enum face,
 axis_enum axis) const throw()
{
  ERROR("BoundaryValue::enforce()",
	"Not implemented yet");
}

//----------------------------------------------------------------------

