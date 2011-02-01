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

  unit_class ("Error");

  //----------------------------------------------------------------------
  printf ("Warning message:\n");

  char warning_message[ERROR_MESSAGE_LENGTH];
  sprintf (warning_message,"Warning message test");
  WARNING_MESSAGE("main",warning_message);

  unit_func("WARNING_MESSAGE");
  unit_assert (true);

  //----------------------------------------------------------------------
  printf ("Incomplete message:\n");

  char incomplete_message[ERROR_MESSAGE_LENGTH];
  sprintf (incomplete_message,"Incomplete message test");
  INCOMPLETE_MESSAGE("main",incomplete_message);

  unit_func("INCOMPLETE_MESSAGE");
  unit_assert (true);

  unit_finalize();

  PARALLEL_EXIT;

}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Error.def.h)
