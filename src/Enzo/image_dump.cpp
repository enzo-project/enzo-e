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

  // slice
  sprintf (filename,"slice-%s-%06d.png",file_root,cycle);

  Monitor * monitor = Monitor::instance();

  if (nz == 1) {
    // 2D: "reduce" along z
    monitor->image(filename,
		   nx,ny,
		   BaryonField[field_density],
		   nx,ny,nz,
		   nx,ny,nz,
		   0,0,0,
		   //		3,3,0,nx-3,ny-3,1,
		   axis_z,reduce_sum, lower/nx, upper/nx);
  } else {
    // 3D projection
    sprintf (filename,"project-%s-%06d-x.png",file_root,cycle);
    monitor->image(filename,
		   ny,nz,
		   BaryonField[field_density],
		   nx,ny,nz,
		   nx,ny,nz,
		   0,0,0,
		   //		  3,3,3,nx-3,ny-3,nz-3,
		   axis_x,reduce_sum,lower, upper);
    sprintf (filename,"project-%s-%06d-y.png",file_root,cycle);
    monitor->image(filename,
		   nz,nx,
		   BaryonField[field_density],
		   nx,ny,nz,
		   nx,ny,nz,
		   0,0,0,
		   //		  3,3,3,nx-3,ny-3,nz-3,
		   axis_y,reduce_sum,lower, upper);
    sprintf (filename,"project-%s-%06d-z.png",file_root,cycle);
    monitor->image(filename,
		   nx,ny,
		   BaryonField[field_density],
		   nx,ny,nz,
		   nx,ny,nz,
		   0,0,0,
		   //		  3,3,3,nx-3,ny-3,nz-3,
		   axis_z,reduce_sum,lower, upper);
  }

}
