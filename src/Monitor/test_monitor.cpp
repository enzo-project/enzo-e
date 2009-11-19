//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      test_monitor.cpp
 * @brief     Program implementing unit tests for monitor classes
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

#include "test.hpp"
#include "monitor.hpp"

int main(int argc, char ** argv)
{
  unit_class ("Monitor");

  unit_open();

  unit_func("Monitor");

  int n = 128;

  // Allocate array

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
	double r = x*x + y*y + z*z;
	if (n/8 < r && r < n/4) array[i] = 1;
      }
    }
  }

  Monitor monitor;

  int map[] = {0,0,0,1,1,1};
  
  monitor.image("test1.png",array,n,n,n,0,0,0,n,n,n,0,reduce_sum,0,1,map,2);

  unit_assert(0);
  unit_close();
}
