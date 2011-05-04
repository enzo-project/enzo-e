// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Parse.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-06
/// @brief    Test program for reading in parameters then displaying them

#include "test.hpp"

#include "parameters.hpp"

#include PARALLEL_CHARM_INCLUDE(test.decl.h)

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;
  cello_parameters_read(stdin);
  cello_parameters_print();
  PARALLEL_EXIT;
}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test.def.h)

