// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Prolong.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implementation of the base Prolong class

#include "problem.hpp"

//----------------------------------------------------------------------

Prolong::Prolong() throw ()
  : monotonic_(false),
    positive_(false)
{
  TRACE("Prolong::Prolong");
}

//======================================================================

