// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      test_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Apr 23 12:25:18 PDT 2008
/// @brief     Program implementing unit tests for performance classes
 
#include "test.hpp"

#include "performance.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Performance.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

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

  unit_func("Timer","start");
  unit_assert((timer.value() - 1.0) < time_tolerance);

  timer.start();
  err = system("sleep 1");
  if (err == -1) ERROR("main","system(sleep) failed!!");
  timer.stop();

  PARALLEL_PRINTF ("Initial timer value = %24.16f\n",timer.value());

  unit_func("Timer","stop");
  unit_assert((timer.value() - 2.0) < time_tolerance);

  unit_func("Performance","Performance");

  Performance * performance = new Performance ();

  // Add attributes

  bool is_monotonic;
  int attribute_timestep =
    performance->new_attribute("timestep", is_monotonic=true);
  int attribute_level = 
    performance->new_attribute("level", is_monotonic=false);
  
  // Add counters

  int counter_time_real = 
    performance->new_counter("time_real");
  int counter_time_sim = 
    performance->new_counter("time_sim");
  int counter_mem_curr_bytes = 
    performance->new_counter("mem_curr_bytes");
  int counter_mem_high_bytes = 
    performance->new_counter("mem_high_bytes");
  int counter_mem_new_count = 
    performance->new_counter("mem_new_count");
  int counter_mem_delete_count = 
    performance->new_counter("mem_delete_count");
  int counter_mem_new_bytes = 
    performance->new_counter("mem_new_bytes");
  int counter_mem_delete_bytes = 
    performance->new_counter("mem_delete_bytes");
  int counter_disk_read_bytes = 
    performance->new_counter("disk_read_bytes");
  int counter_disk_write_bytes = 
    performance->new_counter("disk_write_bytes");
  int counter_disk_read_time = 
    performance->new_counter("disk_read_time");
  int counter_disk_write_time = 
    performance->new_counter("disk_write_time");
  int counter_user_patch_count = 
    performance->new_counter("user_patch_count");
  int counter_user_cell_count = 
    performance->new_counter("user_cell_count");
  int counter_user_particle_count = 
    performance->new_counter("user_particle_count");
  int counter_comm_send_bytes = 
    performance->new_counter("comm_send_bytes");
  int counter_comm_recv_bytes = 
    performance->new_counter("comm_recv_bytes");
  int counter_comm_send_time = 
    performance->new_counter("comm_send_time");
  int counter_comm_recv_time = 
    performance->new_counter("comm_recv_time");
  int counter_comm_global_time = 
    performance->new_counter("comm_global_time");
  int counter_comm_send_count = 
    performance->new_counter("comm_send_count");
  int counter_comm_recv_count = 
    performance->new_counter("comm_recv_count");
  int counter_comm_global_count = 
    performance->new_counter("comm_global_count");

  // Add groups

  int group_1 = performance->new_group("Group 1");
  int group_2 = performance->new_group("Group 2");
  int group_3 = performance->new_group("Group 3");
  
  // Add functions

  int region_1 = performance->new_region("function_1");
  int region_2 = performance->new_region("function_2");
  int region_3 = performance->new_region("function_3");

  // Initialize counters that are non-zero at start

  unit_assert(true);

  unit_func("Performance","new_attribute");
  unit_assert (unit_incomplete);

  unit_func("Performance","new_counter");
  unit_assert (unit_incomplete);

  unit_func("Performance","new_group");
  unit_assert (unit_incomplete);

  unit_func("Performance","new_region");
  unit_assert (unit_incomplete);


  //--------------------------------------------------
  // Attributes
  //--------------------------------------------------

  unit_func("Performance","attribute");
  unit_assert (unit_incomplete);

  unit_func("Performance","set_attribute");
  unit_assert (unit_incomplete);

  unit_func("Performance","num_attributes");
  unit_assert (unit_incomplete);

  //--------------------------------------------------
  // Groups
  //--------------------------------------------------

  unit_func("Performance","begin_group");


  unit_func("Performance","group");
  unit_assert (unit_incomplete);

  unit_func("Performance","set_group");
  unit_assert (unit_incomplete);

  unit_func("Performance","num_groups");
  unit_assert (unit_incomplete);

  unit_func("Performance","end_group");
  unit_assert (unit_incomplete);

  //--------------------------------------------------
  // Regions
  //--------------------------------------------------


  unit_func("Performance","region");
  unit_assert (unit_incomplete);

  unit_func("Performance","set_region");
  unit_assert (unit_incomplete);

  unit_func("Performance","num_regions");
  unit_assert (unit_incomplete);

  unit_func("Performance","start_region");
  unit_assert (unit_incomplete);

  unit_func("Performance","stop_region");
  unit_assert (unit_incomplete);

  //--------------------------------------------------
  // Counters
  //--------------------------------------------------

  unit_func("Performance","counter");
  unit_assert (unit_incomplete);

  unit_func("Performance","set_counter");
  unit_assert (unit_incomplete);

  unit_func("Performance","increment_counter");
  unit_assert (unit_incomplete);

  unit_func("Performance","num_counters");
  unit_assert (unit_incomplete);

  //--------------------------------------------------
  // Disk
  //--------------------------------------------------

  unit_func("Performance","flush");
  unit_assert (unit_incomplete);

  unit_func("Performance","~Performance");
  unit_assert (unit_incomplete);

  delete performance;

  unit_finalize();

  PARALLEL_EXIT;
}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Performance.def.h)
