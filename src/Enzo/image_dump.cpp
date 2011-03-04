// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      image_dump.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Nov 18 2009
/// @brief     Project fields to images

#include "cello.hpp"

#include "enzo.hpp"

void EnzoBlock::image_dump
(
 const char * file_root, 
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

  // slice
  sprintf (filename,"slice-%s-%06d.png",file_root,cycle);

  Monitor * monitor = Monitor::instance();

  monitor->image(filename,
		 BaryonField[field_density],nx,ny,nz,
		 //		3,3,0,nx-3,ny-3,1,
		2,reduce_sum, lower/nx, upper/nx, map,2);

  if (nz > 1) {
    // projection
    sprintf (filename,"project-%s-%06d-x.png",file_root,cycle);
    monitor->image(filename,
		  BaryonField[field_density],nx,ny,nz,
		   //		  3,3,3,nx-3,ny-3,nz-3,
		  0,reduce_sum,lower, upper, map,2);
    sprintf (filename,"project-%s-%06d-y.png",file_root,cycle);
    monitor->image(filename,
		  BaryonField[field_density],nx,ny,nz,
		   //		  3,3,3,nx-3,ny-3,nz-3,
		  1,reduce_sum,lower, upper, map,2);
    sprintf (filename,"project-%s-%06d-z.png",file_root,cycle);
    monitor->image(filename,
		  BaryonField[field_density],nx,ny,nz,
		   //		  3,3,3,nx-3,ny-3,nz-3,
		  2,reduce_sum,lower, upper, map,2);
  }

}
