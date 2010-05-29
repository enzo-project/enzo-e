// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file      image_dump.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Nov 18 2009
/// @brief     Project fields to images

#include "cello.hpp"
#include "cello_hydro.h"

#include <string>

#include "monitor.hpp"

void image_dump(const char * file_root, 
		int cycle, 
		double lower, 
		double upper)
{ 

  int nx = GridDimension[0];
  int ny = GridDimension[1];
  int nz = GridDimension[2];

  char filename[80];

  // color map
  double map[] = {1,1,1, 0,0,0};

  Monitor * monitor = Monitor::instance();

  // slice
  sprintf (filename,"slice-%s-%06d-z.png",file_root,cycle);

  monitor->image(filename,
		BaryonField[field_density],nx,ny,nz,
		3,3,0,nx-3,ny-3,1,
		2,reduce_sum, lower/nx, upper/nx, map,2);

  double min,max;
  min = max = BaryonField[field_density][3 + nx*3];
  for (int ix=3; ix<nx-3; ix++) {
  for (int iy=3; iy<ny-3; iy++) {
    int i = ix + nx*iy;
    double value = BaryonField[field_density][i];
    min = MIN(min,value);
    max = MAX(max,value);
  }
  }
  printf ("Min = %g  Max = %g\n",min,max);
    

  if (nz > 1) {
    // projection
    sprintf (filename,"project-%s-%06d-z.png",file_root,cycle);
    monitor->image(filename,
		  BaryonField[field_density],nx,ny,nz,
		  3,3,3,nx-3,ny-3,nz-3,
		  2,reduce_sum,lower, upper, map,2);
  }

}
