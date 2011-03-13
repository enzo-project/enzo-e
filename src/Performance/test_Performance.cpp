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

  unit_class ("Timer");

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

  unit_class ("Performance");

  unit_func("Performance");

  Performance * performance = new Performance
	  ( num_attributes, num_counters, num_components, 3 );

  // Add attributes

  performance->new_attribute(attribute_timestep, "timestep", 
			     attribute_type_monotonic);
  performance->new_attribute(attribute_level,    "level");
  
  // Add counters

  performance->new_counter(counter_time_real,          "time_real");
  performance->new_counter(counter_time_sim,           "time_sim");
  performance->new_counter(counter_mem_curr_bytes,     "mem_curr_bytes");
  performance->new_counter(counter_mem_high_bytes,     "mem_high_bytes");
  performance->new_counter(counter_mem_new_count,      "mem_new_count");
  performance->new_counter(counter_mem_delete_count,   "mem_delete_count");
  performance->new_counter(counter_mem_new_bytes,      "mem_new_bytes");
  performance->new_counter(counter_mem_delete_bytes,   "mem_delete_bytes");
  performance->new_counter(counter_disk_read_bytes,    "disk_read_bytes");
  performance->new_counter(counter_disk_write_bytes,   "disk_write_bytes");
  performance->new_counter(counter_disk_read_time,     "disk_read_time");
  performance->new_counter(counter_disk_write_time,    "disk_write_time");
  performance->new_counter(counter_user_patch_count,   "user_patch_count");
  performance->new_counter(counter_user_cell_count,    "user_cell_count");
  performance->new_counter(counter_user_particle_count,"user_particle_count");
  performance->new_counter(counter_comm_send_bytes,    "comm_send_bytes");
  performance->new_counter(counter_comm_recv_bytes,    "comm_recv_bytes");
  performance->new_counter(counter_comm_send_time,     "comm_send_time");
  performance->new_counter(counter_comm_recv_time,     "comm_recv_time");
  performance->new_counter(counter_comm_global_time,   "comm_global_time");
  performance->new_counter(counter_comm_send_count,    "comm_send_count");
  performance->new_counter(counter_comm_recv_count,    "comm_recv_count");
  performance->new_counter(counter_comm_global_count,  "comm_global_count");

  // Add groups

  for (int component = component_first;
       component <= num_components;
       component++) {
    performance->new_group(component, component_name[component]);
  }
  
  // Add functions

  performance->new_region(1, "function_1");
  performance->new_region(2, "function_2");
  performance->new_region(3, "function_3");

  // Initialize counters that are non-zero at start

  unit_assert(true);

  unit_func("new_attribute");
  unit_assert (unit_incomplete); //FAILS

  unit_func("new_counter");
  unit_assert (unit_incomplete); //FAILS

  unit_func("new_group");
  unit_assert (unit_incomplete); //FAILS

  unit_func("new_region");
  unit_assert (unit_incomplete); //FAILS


  //--------------------------------------------------
  // Attributes
  //--------------------------------------------------

  unit_func("attribute");
  unit_assert (unit_incomplete); //FAILS

  unit_func("set_attribute");
  unit_assert (unit_incomplete); //FAILS

  unit_func("num_attributes");
  unit_assert (unit_incomplete); //FAILS

  //--------------------------------------------------
  // Groups
  //--------------------------------------------------

  unit_func ("begin_group");


  unit_func("group");
  unit_assert (unit_incomplete); //FAILS

  unit_func("set_group");
  unit_assert (unit_incomplete); //FAILS

  unit_func("num_groups");
  unit_assert (unit_incomplete); //FAILS

  unit_func("end_group");
  unit_assert (unit_incomplete); //FAILS

  //--------------------------------------------------
  // Regions
  //--------------------------------------------------


  unit_func("region");
  unit_assert (unit_incomplete); //FAILS

  unit_func("set_region");
  unit_assert (unit_incomplete); //FAILS

  unit_func("num_regions");
  unit_assert (unit_incomplete); //FAILS

  unit_func("start_region");
  unit_assert (unit_incomplete); //FAILS

  unit_func("stop_region");
  unit_assert (unit_incomplete); //FAILS

  //--------------------------------------------------
  // Counters
  //--------------------------------------------------

  unit_func("counter");
  unit_assert (unit_incomplete); //FAILS

  unit_func("set_counter");
  unit_assert (unit_incomplete); //FAILS

  unit_func("increment_counter");
  unit_assert (unit_incomplete); //FAILS

  unit_func("num_counters");
  unit_assert (unit_incomplete); //FAILS

  //--------------------------------------------------
  // Disk
  //--------------------------------------------------

  unit_func("flush");
  unit_assert (unit_incomplete); //FAILS

  unit_func("~Performance");
  unit_assert (unit_incomplete); //FAILS

  delete performance;


  Papi papi;

  papi.start();
  err = system("sleep 1");
  if (err == -1) ERROR("main","system(sleep) failed!!");
  papi.stop();

  PARALLEL_PRINTF ("time_real = %f\n",papi.time_real());
  PARALLEL_PRINTF ("time_proc = %f\n",papi.time_proc());

  papi.start();
  float a=1.0,b=2.5,c=0;
  c = a*b+3.5;
  papi.stop();

  PARALLEL_PRINTF ("c=%f\n",c);
  PARALLEL_PRINTF ("time_real  = %f\n",papi.time_real());
  PARALLEL_PRINTF ("time_proc  = %f\n",papi.time_proc());
  PARALLEL_PRINTF ("flop_count = %lld\n",papi.flop_count());
  PARALLEL_PRINTF ("flop_rate = %f\n",papi.flop_rate());


  unit_finalize();

  PARALLEL_EXIT;
}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Performance.def.h)
