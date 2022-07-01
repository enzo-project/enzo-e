// See LICENSE_CELLO file for license and copyright information

/// @file     test_Parse.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-06
/// @brief    Test program for reading in parameters then displaying them

#include "main.hpp" 
#include "test.hpp"

#include "parameters.hpp"

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;
  unit_init (0, 1);
  unit_class("Parameters");
  if (PARALLEL_ARGC < 2){
    CkPrintf("ERROR: This test expects a parameter file to be passed as a "
             "command line argument. None were provided.\n");
    unit_assert(false);
  } else {
    for (int i=1; i<PARALLEL_ARGC; i++) {
      const char * filename = PARALLEL_ARGV[i];
      unit_func("parse");
      FILE * fp = fopen (filename,"r");
      cello_parameters_read(filename,fp);
      unit_assert(true);
    }
    cello_parameters_print();
  }

  unit_finalize();
  exit_();
}
PARALLEL_MAIN_END
