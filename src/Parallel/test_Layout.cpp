// $Id: test_block.cpp 1369 2010-04-08 01:38:06Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-19
/// @brief    Unit tests for the Layout class
///
/// Run with mpirun -np 4


#include <mpi.h>

#include "cello.h"

#include "error.hpp"
#include "test.hpp"
#include "parallel.hpp"

#define INDEX3(I,N) I[0] + N[0]*(I[1] + N[1]*I[2])

int main(int argc, char ** argv)
{
  unit_class("Layout");
  unit_func("Layout");
  Layout layout_serial;
  unit_assert (true);

  layout_serial.set_periodic(axis_x,true);
  layout_serial.set_periodic(axis_y,true);
  layout_serial.set_periodic(axis_z,true);

  unit_func("process_count");
  unit_assert (layout_serial.process_count() == 1);
  unit_func("thread_count");
  unit_assert (layout_serial.thread_count()  == 1);
  unit_func("data_blocks_per_process");
  unit_assert (layout_serial.data_blocks_per_process()  == 1);
  unit_func("data_blocks_per_thread");
  unit_assert (layout_serial.data_blocks_per_thread()  == 1);
  unit_func("is_periodic");
  unit_assert (layout_serial.is_periodic(axis_x) == true);
  unit_assert (layout_serial.is_periodic(axis_y) == true);
  unit_assert (layout_serial.is_periodic(axis_z) == true);

//   unit_assert (neighbor_is_internal(0,0,0,axis_x,face_upper);
//   unit_assert (neighbor_is_internal(0,0,0,axis_x,face_lower);
//   unit_assert (neighbor_is_internal(0,0,0,axis_y,face_upper);
//   unit_assert (neighbor_is_internal(0,0,0,axis_y,face_lower);
//   unit_assert (neighbor_is_internal(0,0,0,axis_z,face_upper);
//   unit_assert (neighbor_is_internal(0,0,0,axis_z,face_lower);

  int nx=5, ny=3, nz=1;
  bool passed = true;
  for (int iz=0; iz<nz && passed; iz++) {
    for (int iy=0; iy<ny && passed; iy++) {
      for (int ix=0; ix<nx && passed; ix++) {

	int i,j,jx,jy,jz;

	index_3_to_1(j,ix,iy,iz,nx,ny,nz);
	//	jz=(((((j-j%nx)/nx)-(jy=((j-jx=(j%nx))/nx)%ny))/ny)%nz);
	index_1_to_3(j,jx,jy,jz,nx,ny,nz);
	passed = passed && (ix==jx && iy==jy && iz==jz);

	i = ix + nx*(iy + ny*iz);

	index_1_to_3(i,jx,jy,jz,nx,ny,nz);
	index_3_to_1(j,jx,jy,jz,nx,ny,nz);
	passed = passed && (i == j);
	
      }
    }
  }
  unit_assert (passed);
}
