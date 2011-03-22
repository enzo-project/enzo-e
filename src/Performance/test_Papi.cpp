// $Id: test_Papi.cpp 1696 2010-08-04 05:56:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_TEMPLATE.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Papi class

#include "test.hpp"

#include "performance.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Papi.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class("Papi");

  Papi papi;

  papi.start();
  int err = system("sleep 1");
  if (err == -1) ERROR("main","system(sleep) failed!!");
  float a=1.0,b=2.5,c=0;
  c = a*b+3.5+sin(3.0);
  papi.stop();

  PARALLEL_PRINTF ("c=%f\n",c); // prevent compiler optimizing computations out
  unit_func("time_real");
  unit_assert(papi.time_real() > 0);
  unit_func("time_proc");
  unit_assert(papi.time_proc() > 0);
  unit_func("flop_count");
  unit_assert(papi.flop_count() > 0);
  unit_func("flop_rate");
  unit_assert(papi.flop_rate() > 0);
  
  PARALLEL_PRINTF ("time_real  = %f\n",papi.time_real());
  PARALLEL_PRINTF ("time_proc  = %f\n",papi.time_proc());
  PARALLEL_PRINTF ("flop_count = %lld\n",papi.flop_count());
  PARALLEL_PRINTF ("flop_rate = %f\n",papi.flop_rate());

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Papi.def.h)
