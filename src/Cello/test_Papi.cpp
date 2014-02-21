// See LICENSE_CELLO file for license and copyright information

/// @file     test_TEMPLATE.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Papi class

#include "main.hpp" 
#include "test.hpp"

#include "performance.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Papi");

  Papi papi;

  papi.init();

  //  papi.add_event(PAPI_FP_INS);
  papi.add_event("PAPI_FP_OPS");

  const int num_events = papi.num_events();

  const int num_count = 4;

  // - count means turn off PAPI counting: result should be 0
  const int count_array[] = {10000,1000,-1000,1000};

  papi.start_events();

  long long * values_start = new long long [num_events];
  long long * values_stop  = new long long [num_events];

  for (int index_count = 0; index_count<num_count; index_count++) {

    int count = abs(count_array[index_count]);

    if (count_array[index_count] < 0) {
      papi.stop_events();
    }

    papi.event_values(values_start);

    float a=1.0, b=2.5;
    for (int i=0; i<count; i++) {
      b = a + b;
    }

    papi.event_values(values_stop);

    if (count_array[index_count] < 0) {
      papi.start_events();
    }

    for (int ie = 0; ie < num_events; ie++) {
      values_stop[ie] -= values_start[ie];
      printf ("event %s value %lld  [inhibit optimize: b = %f]\n",
	      papi.event_name(ie).c_str(),
	      values_stop[ie],b);
    }

    if (count_array[index_count] > 0) 
      unit_assert (fabs(values_stop[0] - count)/(count) < 0.05);
    else
      unit_assert (values_stop[0] == 0);

  }

  papi.stop_events();

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
