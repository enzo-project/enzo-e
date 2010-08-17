// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      test_Schedule.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:47:35 PST 2008
/// @brief     Program implementing unit tests for the Schedule class
 
#include "cello.hpp"

#include "error.hpp"
#include "test.hpp"

#include "parallel.def"

#include PARALLEL_CHARM_INCLUDE(test_Schedule.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  bool passed = false;

  unit_class ("Schedule");
  unit_func ("Schedule");

  //  Schedule schedule;

  //FAILS
  unit_assert(passed);

  unit_finalize();

  PARALLEL_EXIT;

}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Schedule.def.h)
