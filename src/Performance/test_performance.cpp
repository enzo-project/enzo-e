//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      test_performance.cpp
 * @brief     Program implementing unit tests for performance classes
 * @author    James Bordner
 * @date      Wed Apr 23 12:25:18 PDT 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "error.hpp"
#include "test.hpp"
#include "array.hpp"
#include "performance.hpp"

int main(int argc, char ** argv)
{
  Timer timer;
  const double time_tolerance = 0.05;

  // Timer tests

  unit_class ("Timer");

  unit_open();

  printf ("Initial timer value = %24.16f\n",timer.value());

  int return_value;
  timer.start();

  return_value = system("sleep 1");
  timer.stop();

  printf ("Initial timer value = %24.16f\n",timer.value());

  unit_assert((timer.value() - 1.0) < time_tolerance);

  timer.start();
  return_value = system("sleep 1");
  timer.stop();

  printf ("Initial timer value = %24.16f\n",timer.value());

  unit_assert((timer.value() - 2.0) < time_tolerance);

  unit_class ("Performance");

  unit_func("Performance");
  Performance performance ( num_perf_attributes,
			    num_perf_counters,
			    num_components,
			    3 ); // 3 regions
  unit_assert(true);

  //--------------------------------------------------
  // Attributes
  //--------------------------------------------------

  unit_func("new_attribute");
  unit_assert (false);

  unit_func("get_attribute");
  unit_assert (false);

  unit_func("set_attribute");
  unit_assert (false);

  unit_func("num_attributes");
  unit_assert (false);

  //--------------------------------------------------
  // Groups
  //--------------------------------------------------

  unit_func ("begin_group");


  unit_func("new_group");
  unit_assert (false);

  unit_func("get_group");
  unit_assert (false);

  unit_func("set_group");
  unit_assert (false);

  unit_func("num_groups");
  unit_assert (false);

  unit_func("end_group");
  unit_assert (false);

  //--------------------------------------------------
  // Regions
  //--------------------------------------------------


  unit_func("new_region");
  unit_assert (false);

  unit_func("get_region");
  unit_assert (false);

  unit_func("set_region");
  unit_assert (false);

  unit_func("num_regions");
  unit_assert (false);

  unit_func("start_region");
  unit_assert (false);

  unit_func("stop_region");
  unit_assert (false);

  //--------------------------------------------------
  // Counters
  //--------------------------------------------------

  unit_func("new_counter");
  unit_assert (false);

  unit_func("get_counter");
  unit_assert (false);

  unit_func("set_counter");
  unit_assert (false);

  unit_func("increment_counter");
  unit_assert (false);

  unit_func("num_counters");
  unit_assert (false);

  //--------------------------------------------------
  // Disk
  //--------------------------------------------------

  unit_func("flush");
  unit_assert (false);

  unit_func("~Performance");
  unit_assert (false);

  unit_close();
}
