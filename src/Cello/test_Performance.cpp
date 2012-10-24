// See LICENSE_CELLO file for license and copyright information

/// @file      test_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Apr 23 12:25:18 PDT 2008
/// @brief     Program implementing unit tests for performance classes

#include "main.hpp" 
#include "test.hpp"

#include "performance.hpp"

//----------------------------------------------------------------------

void sleep_flop (int s, int count)
{
  char sleep_string [10];
  sprintf (sleep_string,"sleep %d",s);
  int err = system(sleep_string);

  if (err == -1) ERROR("main","system(sleep) failed!!");

  float a=1.0, b=2.5;
  for (int i=0; i<count; i++) {
    b = a + b;
  }
  printf ("inhibit optimizing out loop %f\n",b);

}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Performance");

  Performance * performance = new Performance ();

  // Initialize counters that are non-zero at start

  unit_assert(performance != NULL);

  unit_func("new_counter");

  int index_counter_1 = performance->new_counter("counter_1");
  int index_counter_2 = performance->new_counter("counter_2");

  unit_assert (index_counter_1 != index_counter_2);

  unit_func("num_counters");

  int num_counters = performance->num_counters();
  unit_assert (num_counters >= 2);

  printf ("num counters = %d\n",num_counters);
  long long * counters = new long long [num_counters];

  unit_func("new_region");

  int index_region_1 = performance->new_region("region_1");
  int index_region_2 = performance->new_region("region_2");

  unit_assert (index_region_1 != index_region_2);

  performance->begin();

  performance->start_region(index_region_1);

  sleep_flop (1,1000000);
  performance->increment_counter(index_counter_1,10);
  performance->increment_counter(index_counter_2,20);

  performance->start_region(index_region_2);
  
  sleep_flop (2,500000);
  performance->increment_counter(index_counter_1,50);
  performance->increment_counter(index_counter_2,100);

  performance->stop_region(index_region_2);

  sleep_flop (1,1000000);
  performance->increment_counter(index_counter_1,10);
  performance->increment_counter(index_counter_2,20);

  performance->stop_region(index_region_1);

  unit_func("user counters");

  performance->region_counters(index_region_1,counters);

  printf ("%lld %lld %lld\n",counters[0],counters[1],counters[2]);
  unit_assert(counters[index_counter_1] == 10 + 50 + 10);
  unit_assert(counters[index_counter_2] == 20 + 100 + 20);

  performance->region_counters(index_region_2,counters);
  printf ("%lld %lld %lld\n",counters[0],counters[1],counters[2]);

  unit_assert(counters[index_counter_1] == 50);
  unit_assert(counters[index_counter_2] == 100);

  performance->end();

  int num_regions = performance->num_regions();
  for (int ir = 0; ir < num_regions; ir++) {

    performance->region_counters(ir,counters);

    for (int ic = 0; ic < num_counters; ic++) {
    
      printf ("region %s counter %s value %lld\n",
	      performance->region_name(ir).c_str(),
	      performance->counter_name(ic).c_str(),
	      counters[ic]);
    }
  }

  delete performance;

  unit_finalize();

  PARALLEL_EXIT;
}
PARALLEL_MAIN_END
