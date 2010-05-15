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

  //----------------------------------------------------------------------
  // index conversions
  //----------------------------------------------------------------------

  unit_func("index_*");

  int nx=5, ny=3, nz=4;
  bool passed = true;
  for (int iz=0; iz<nz && passed; iz++) {
    for (int iy=0; iy<ny && passed; iy++) {
      for (int ix=0; ix<nx && passed; ix++) {

	int i,j,jx,jy,jz,kx,ky,kz;

	// ix,iy,iz -> i
	j = index_3_to_1(ix,iy,iz,nx,ny,nz);
	index_1_to_3(j,jx,jy,jz,nx,ny,nz);
	passed = passed && (ix==jx && iy==jy && iz==jz);

	// i -> ix, i -> iy, i -> iz
	kx = index_1_to_x(j,nx,ny,nz);
	ky = index_1_to_y(j,nx,ny,nz);
	kz = index_1_to_z(j,nx,ny,nz);
	passed = passed && (ix==kx && iy==ky && iz==kz);

	i = ix + nx*(iy + ny*iz);

	// i -> ix,iy,iz
	index_1_to_3(i,jx,jy,jz,nx,ny,nz);
	j = index_3_to_1(jx,jy,jz,nx,ny,nz);
	passed = passed && (i == j);
	
      }
    }
  }
  unit_assert (passed);

  //----------------------------------------------------------------------
  // serial layout: (processes,threads,data blocks) = (1,1,1)
  //----------------------------------------------------------------------

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

  // periodic, so all neighbors should be internal
  unit_assert (layout_serial.neighbor_is_internal(0,0,0,axis_x,+1));
  unit_assert (layout_serial.neighbor_is_internal(0,0,0,axis_x,-1));
  unit_assert (layout_serial.neighbor_is_internal(0,0,0,axis_y,+1));
  unit_assert (layout_serial.neighbor_is_internal(0,0,0,axis_y,-1));
  unit_assert (layout_serial.neighbor_is_internal(0,0,0,axis_z,+1));
  unit_assert (layout_serial.neighbor_is_internal(0,0,0,axis_z,-1));

}
