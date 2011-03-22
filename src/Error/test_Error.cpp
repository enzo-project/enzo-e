// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      test_Error.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Aug 20 11:24:14 PDT 2008
/// @brief     Program implementing unit tests for error classes
 
#include "test.hpp"

#include "error.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Error.decl.h)

PARALLEL_MAIN_BEGIN

{

  PARALLEL_INIT;

  unit_init();

  //----------------------------------------------------------------------
  PARALLEL_PRINTF ("Warning message:\n");

  char warning_message[ERROR_LENGTH];
  sprintf (warning_message,"Warning message test");
  WARNING("main",warning_message);

  unit_func("Error","WARNING");
  unit_assert (true);

  //----------------------------------------------------------------------
  PARALLEL_PRINTF ("Incomplete message:\n");

  INCOMPLETE("main");

  unit_func("Error","INCOMPLETE");
  unit_assert (true);

  unit_finalize();

  PARALLEL_EXIT;

}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Error.def.h)
