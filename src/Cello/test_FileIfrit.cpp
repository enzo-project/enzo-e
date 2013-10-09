// See LICENSE_CELLO file for license and copyright information

/// @file     test_FileIfrit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:47:35 PST 2008
/// @brief    Program implementing unit tests for the FileIfrit class
 
#include "main.hpp" 
#include "test.hpp"

#include "disk.hpp"
PARALLEL_MAIN_BEGIN

{
  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("FileIfrit");

  int nx = 64;
  int ny = 64;
  int nz = 64;

  const char filename[] = "test_disk.bin";
  int n = nx*ny*nz;

  float * a = new float[n];

  for (int iz=0; iz<nz; iz++) {
    double z = 1.0*iz/(nz-1);
    for (int iy=0; iy<ny; iy++) {
      double y = 1.0*iy/(ny-1);
      for (int ix=0; ix<nx; ix++) {
	double x = 1.0*ix/(nx-1);
	int i = ix + nx*(iy + ny*iz);
	a[i] = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5);
      }
    }
  }

  FileIfrit ifrit;

  unit_func("write_bin");
  ifrit.write_bin(filename,a,nx,ny,nz);
  unit_assert(true);

  unit_func("read_bin");
  float * b = new float[n];
  int mx,my,mz;
  ifrit.read_bin(filename,b,&mx,&my,&mz);
  unit_assert (mx == nx);
  unit_assert (my == ny);
  unit_assert (mz == nz);

  bool passed = true;
  for (int iz=0; iz<nz && passed; iz++) {
    for (int iy=0; iy<ny && passed; iy++) {
      for (int ix=0; ix<nx && passed; ix++) {
	int i = ix + nx*(iy + ny*iz);
	if (a[i] != b[i]) passed = false;
      }
    }
  }

  delete [] a;
  delete [] b;

  unit_assert(passed);

  unit_finalize();

  PARALLEL_EXIT;
}
PARALLEL_MAIN_END
