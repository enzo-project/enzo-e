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

  Parameters parameters;
  unit_class("Performance");

  Performance * performance = new Performance (NULL);

  // Initialize counters that are non-zero at start

  unit_assert(performance != NULL);

  unit_func("new_counter");

  int id_counter_1 = performance->new_counter(counter_type_user,"counter_1");
  int id_counter_2 = performance->new_counter(counter_type_user,"counter_2");
  int id_counter_flops =
    performance->new_counter(counter_type_papi,"PAPI_FP_INS");

  unit_assert (id_counter_1 != id_counter_2);
  unit_assert (id_counter_1 != id_counter_flops);
  unit_assert (id_counter_2 != id_counter_flops);

  unit_func("num_counters");

  int num_counters = performance->num_counters();
  unit_assert (num_counters >= 2);

  printf ("num counters = %d\n",num_counters);
  long long * region_counters = new long long [num_counters];

  unit_func("new_region");

  int id_region_1 = 0;
  int id_region_2 = 1;
  performance->new_region(id_region_1,"region_1");
  performance->new_region(id_region_2,"region_2");

  unit_assert (id_region_1 != id_region_2);

  performance->begin();

  long long counter_values [10];
  TRACE2("%d %d",performance->counter_values(counter_values) , num_counters);
  unit_assert (performance->counter_values(counter_values) == num_counters);

  performance->start_region(id_region_1);


  sleep_flop (1,1000000);
  performance->increment_counter(id_counter_1,10);

  performance->increment_counter(id_counter_2,20);


  performance->start_region(id_region_2);

  
  sleep_flop (2,500000);

  performance->increment_counter(id_counter_1,50);

  performance->stop_region(id_region_2);
  performance->start_region(id_region_2);

  performance->increment_counter(id_counter_2,100);


  performance->stop_region(id_region_2);


  sleep_flop (1,1000000);
  performance->increment_counter(id_counter_1,10);

  performance->increment_counter(id_counter_2,20);



  performance->stop_region(id_region_1);

  unit_func("user counters");



  performance->region_counters(id_region_1,region_counters);

  int index_counter_1 = id_counter_1;
  int index_counter_2 = id_counter_2;

  
  unit_assert(region_counters[index_counter_1] == 10 + 50 + 10);
  unit_assert(region_counters[index_counter_2] == 20 + 100 + 20);

  performance->region_counters(id_region_2,region_counters);

  printf ("%lld %d\n",region_counters[index_counter_1] , 50);
  unit_assert(region_counters[index_counter_1] == 50);
  unit_assert(region_counters[index_counter_2] == 100);

  performance->end();

  performance->counter_values(counter_values) ;

  for (int index_counter = 0; index_counter < num_counters; index_counter++) {
    
    int id_counter = index_counter;

    printf ("COUNTER %s VALUE %lld\n",
	    performance->counter_name(id_counter).c_str(),
	    counter_values[index_counter]);
  }

  int num_regions = performance->num_regions();
  for (int ir = 0; ir < num_regions; ir++) {

    performance->region_counters(ir,region_counters);

    for (int ic = 0; ic < num_counters; ic++) {
    
      printf ("region %s value %lld\n",
	      performance->region_name(ir).c_str(),
	      region_counters[ic]);
    }
  }

  delete performance;

  unit_finalize();

  exit_();
}
PARALLEL_MAIN_END
