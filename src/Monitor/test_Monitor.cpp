// See LICENSE_CELLO file for license and copyright information

/// @file     test_Monitor.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-18
/// @brief    Program implementing unit tests for the Monitor class

#include "test.hpp"

#include "monitor.hpp"

#ifdef CONFIG_USE_CHARM
#   include "main.decl.h"
#endif

PARALLEL_MAIN_BEGIN
{
  
  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  Monitor * monitor  = Monitor::instance();

  unit_init(parallel->rank(),parallel->size());

  unit_class("Monitor");

  unit_func("Monitor");
  unit_assert(true);


  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#ifdef CONFIG_USE_CHARM
#   include "main.def.h"
#endif
