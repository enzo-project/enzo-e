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

  //  papi.add_counter(PAPI_FP_INS);
  papi.add_counter("PAPI_FP_OPS");

  int count_array[] = {10000,1000,100};

  for (int icount = 0; icount<2; icount++) {

    int count = count_array[icount];

    papi.start();

    float a=1.0, b=2.5;
    for (int i=0; i<count; i++) {
      b = a + b;
    }
    printf ("printf to inhibit optimizing-out code: flops = %f\n",b);

    papi.stop();

    unit_assert (fabs(papi.value(0) - count)/(count) < 0.05);

    int num_counters = papi.num_counters();
    for (int i=0; i<num_counters; i++) {
      printf ("%d %s %lld\n",i,papi.name(i).c_str(),papi.value(i));
    }
  }



  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
