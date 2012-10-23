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

  papi.add_region("A");
  papi.add_region("B");
  papi.add_region("C");

  int count_array[] = {10000,1000,100};
  std::string region_array[] = {"A","B","C"};

  TRACE0;
  for (int index_region = 0; index_region<3; index_region++) {

    int count = count_array[index_region];

  TRACE0;
    papi.start_region(index_region);
  TRACE0;

    float a=1.0, b=2.5;
    for (int i=0; i<count; i++) {
      b = a + b;
    }
    printf ("printf to inhibit optimizing-out code: flops = %f\n",b);

  TRACE0;
    papi.stop_region(index_region);
  TRACE0;

    const long long * values = papi.values(index_region);
  TRACE0;

  TRACE1 ("values = %p",values);
    unit_assert (fabs(values[0] - count)/(count) < 0.05);

  }

  TRACE0;
  const int num_events = papi.num_events();
  TRACE0;
  const int num_regions = papi.num_regions();
  TRACE0;
  for (int index_event=0; index_event<num_events; index_event++) {

    const long long * events = papi.values(index_event);

    for (int index_region=0; index_region<num_regions; index_region++) {

  TRACE0;
      printf ("region %s  event %s value %lld\n",
	      papi.region_name(index_region).c_str(),
	      papi.event_name(index_event).c_str(),
	      events[index_event]);

    }
  }


  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
