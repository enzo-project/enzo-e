// $Id: test_Papi.cpp 1696 2010-08-04 05:56:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_TEMPLATE.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Papi class

#include "test.hpp"

#include "performance.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Papi.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class ("Papi");

  Papi papi;

  papi.start();
  papi.stop();
 
  unit_assert(papi.time_real() > 0);
  unit_assert(papi.time_proc() > 0);
  unit_assert(papi.flop_count() > 0);
  unit_assert(papi.flop_rate() > 0);
  
  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Papi.def.h)
