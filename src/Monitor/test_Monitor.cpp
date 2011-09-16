// See LICENSE_CELLO file for license and copyright information

/// @file     test_Monitor.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-18
/// @brief    Program implementing unit tests for the Monitor class

#include "main.hpp" 
#include "test.hpp"

#include "monitor.hpp"

PARALLEL_MAIN_BEGIN
{
  
  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  Monitor * monitor  = Monitor::instance();

  unit_init(parallel->rank(),parallel->size());

  unit_class("Monitor");

  unit_func("Monitor");

  unit_assert(monitor != NULL);


  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
