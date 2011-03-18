// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Monitor.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-18
/// @brief    Program implementing unit tests for the Monitor class

#include "test.hpp"

#include "monitor.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Monitor.decl.h)

PARALLEL_MAIN_BEGIN
{
  
  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  Monitor * monitor  = Monitor::instance();

  unit_init(parallel->rank(),parallel->size());

  unit_class ("Monitor");

  int n = 128;

  // Allocate array

  if (parallel->rank()==0) {
    PARALLEL_PRINTF ("pngwriter version = %g\n",pngwriter::version());
  }

  float * array = new float [n*n*n];
  for (int i=0; i<n*n*n; i++) array[i]=0;

  // Set values to 1 for radius between n/8 and n/4 (note n/2 is boundary)
  for (int iz=0; iz<n; iz++) {
    double z = (iz - 0.5*n) / n;
    for (int iy=0; iy<n; iy++) {
      double y = (iy - 0.5*n) / n;
      for (int ix=0; ix<n; ix++) {
  	int i = ix + n*(iy + n*iz);
  	double x = (ix - 0.5*n) / n;
  	//	array[i] = y;
  	double r = sqrt(x*x + y*y + z*z);
  	if (0.25 < r && r < 0.5) array[i] = 1;
      }
    }
  }

  unit_func("Monitor");
  unit_assert(true);

  unit_func("image");

  {
    int axis=0;
    monitor->image("monitor_image_1.png",array,n,n,n,0,0,0, axis,reduce_sum,0,1);
    unit_assert(true);
  }

  {
    double map_r[] = {0.0, 1.0};
    double map_g[] = {0.0, 0.5};
    double map_b[] = {0.5, 1.0};
    int axis=1;
    monitor->set_image_map(2,map_r,map_g,map_b);
    monitor->image("monitor_image_2.png",array,n,n,n,0,0,0, axis,reduce_sum,0,1);
    unit_assert(true);
  }

  {
    double map_r[] = {0.0, 1.0, 0.0, 0.0};
    double map_g[] = {0.0, 0.0, 1.0, 0.0};
    double map_b[] = {0.0, 0.0, 0.0, 1.0};
    int axis=2;
    monitor->set_image_map(4,map_r,map_g,map_b);
    monitor->image("monitor_image_3.png",array,n,n,n,0,0,0, axis,reduce_sum,0,1);
    unit_assert(true);
  }

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Monitor.def.h)
