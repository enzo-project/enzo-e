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
  Performance performance ( num_attributes,
			    num_counters,
			    num_components,
			    3 ); // 3 regions

  // Add attributes

  performance.new_attribute(time_real,          "time_real");
  performance.new_attribute(time_sim,           "time_sim");
  performance.new_attribute(mem_curr_bytes,     "mem_curr_bytes");
  performance.new_attribute(mem_high_bytes,     "mem_high_bytes");
  performance.new_attribute(mem_new_count,      "mem_new_count");
  performance.new_attribute(mem_delete_count,   "mem_delete_count");
  performance.new_attribute(mem_new_bytes,      "mem_new_bytes");
  performance.new_attribute(mem_delete_bytes,   "mem_delete_bytes");
  performance.new_attribute(disk_read_bytes,    "disk_read_bytes");
  performance.new_attribute(disk_write_bytes,   "disk_write_bytes");
  performance.new_attribute(disk_read_time,     "disk_read_time");
  performance.new_attribute(disk_write_time,    "disk_write_time");
  performance.new_attribute(user_patch_count,   "user_patch_count");
  performance.new_attribute(user_cell_count,    "user_cell_count");
  performance.new_attribute(user_particle_count,"user_particle_count");

#ifdef CONFIG_USE_MPI
  performance.new_attribute(comm_send_bytes,    "comm_send_bytes");
  performance.new_attribute(comm_recv_bytes,    "comm_recv_bytes");
  performance.new_attribute(comm_send_time,     "comm_send_time");
  performance.new_attribute(comm_recv_time,     "comm_recv_time");
  performance.new_attribute(comm_global_time,   "comm_global_time");
  performance.new_attribute(comm_send_count,    "comm_send_count");
  performance.new_attribute(comm_recv_count,    "comm_recv_count");
  performance.new_attribute(comm_global_count,  "comm_global_count");
#endif /* CONFIG_USE_MPI */

  // Add counters
  // Add groups
  // Add functions

  // Initialize counters that are non-zero at start

  unit_assert(true);

  unit_func("new_attribute");
  unit_assert (false);

  unit_func("new_counter");
  unit_assert (false);

  unit_func("new_group");
  unit_assert (false);

  unit_func("new_region");
  unit_assert (false);


  //--------------------------------------------------
  // Attributes
  //--------------------------------------------------

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
