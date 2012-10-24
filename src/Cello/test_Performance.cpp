// See LICENSE_CELLO file for license and copyright information

/// @file      test_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Apr 23 12:25:18 PDT 2008
/// @brief     Program implementing unit tests for performance classes

#include "main.hpp" 
#include "test.hpp"

#include "performance.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Timer");

  Timer timer;
  const double time_tolerance = 0.05;

  // Timer tests

  PARALLEL_PRINTF ("Initial timer value = %24.16f\n",timer.value());

  timer.start();

  int err; // To inhibit warnings--I'm not expecting the sleep command to fail

  err = system("sleep 1");
  if (err == -1) ERROR("main","system(sleep) failed!!");
  timer.stop();

  PARALLEL_PRINTF ("Initial timer value = %24.16f\n",timer.value());

  unit_func("start");
  unit_assert((timer.value() - 1.0) < time_tolerance);

  timer.start();
  err = system("sleep 1");
  if (err == -1) ERROR("main","system(sleep) failed!!");
  timer.stop();

  PARALLEL_PRINTF ("Initial timer value = %24.16f\n",timer.value());

  unit_func("stop");
  unit_assert((timer.value() - 2.0) < time_tolerance);

  unit_class("Performance");

  unit_func("Performance");

  Performance * performance = new Performance ();

  // Initialize counters that are non-zero at start

  unit_assert(true);

  unit_func("new_counter");
  unit_assert (unit_incomplete);

  unit_func("new_region");
  unit_assert (unit_incomplete);

  //--------------------------------------------------
  // Regions
  //--------------------------------------------------


  unit_func("region");
  unit_assert (unit_incomplete);

  unit_func("set_region");
  unit_assert (unit_incomplete);

  unit_func("num_regions");
  unit_assert (unit_incomplete);

  unit_func("start_region");
  unit_assert (unit_incomplete);

  unit_func("stop_region");
  unit_assert (unit_incomplete);

  //--------------------------------------------------
  // Counters
  //--------------------------------------------------

  unit_func("counter");
  unit_assert (unit_incomplete);

  unit_func("set_counter");
  unit_assert (unit_incomplete);

  unit_func("increment_counter");
  unit_assert (unit_incomplete);

  unit_func("num_counters");
  unit_assert (unit_incomplete);

  //--------------------------------------------------
  // Disk
  //--------------------------------------------------

  unit_func("flush");
  unit_assert (unit_incomplete);

  unit_func("~Performance");
  unit_assert (unit_incomplete);

  delete performance;

  unit_finalize();

  PARALLEL_EXIT;
}
PARALLEL_MAIN_END
