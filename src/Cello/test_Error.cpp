// See LICENSE_CELLO file for license and copyright information

/// @file      test_Error.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Aug 20 11:24:14 PDT 2008
/// @brief     Program implementing unit tests for error classes
 
#include "main.hpp" 
#include "test.hpp"

#include "error.hpp"

PARALLEL_MAIN_BEGIN

{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Error");
  //----------------------------------------------------------------------
  PARALLEL_PRINTF ("Warning message:\n");

  char warning_message[ERROR_LENGTH];
  sprintf (warning_message,"Warning message test");
  WARNING("main",warning_message);

  unit_func("WARNING");
  unit_assert (true);

  //----------------------------------------------------------------------
  PARALLEL_PRINTF ("Incomplete message:\n");

  INCOMPLETE("main");

  unit_func("INCOMPLETE");
  unit_assert (true);

  unit_finalize();

  exit_();

}

PARALLEL_MAIN_END

