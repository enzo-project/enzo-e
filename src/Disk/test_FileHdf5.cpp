// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_FileHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:47:35 PST 2008
/// @todo     Multiple datasets with different precisions
/// @brief    Program implementing unit tests for the FileHdf5 class
 
#include "test.hpp"

#include "disk.hpp"

// #include "error.hpp"
// #include "test.hpp"
// #include "disk.hpp"

// #include "parallel.def"

#include PARALLEL_CHARM_INCLUDE(test_FileHdf5.decl.h)

PARALLEL_MAIN_BEGIN

{
  PARALLEL_INIT;

  unit_init();

  int nx = 100;
  int ny = 50;
  Scalar * a = new Scalar[nx*ny];

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {
      int i = ix + nx*(iy);
      a[i] = ix*3 + iy*5;
    }
  }

  unit_class ("FileHdf5");

  unit_func("file_open");

  FileHdf5 hdf5;

  int mx,my,mz;
  hdf5.file_open("file_open_test.h5","w");
  mx = nx;
  my = ny;
  mz = 1;

  hdf5.dataset_open_write ("dataset",precision_default,mx,my,mz);
  hdf5.write((char *)a,precision_default);
  hdf5.dataset_close ();
  hdf5.file_close();


  hdf5.file_open("file_open_test.h5","r");
  hdf5.dataset_open_read ("dataset",&mx,&my,&mz);

  Scalar * b = new Scalar[nx*ny];
  
  hdf5.read((char *)b,precision_default);
  hdf5.dataset_close ();
  hdf5.file_close();

  bool passed = true;
  for (int iy=0; iy<ny && passed; iy++) {
    for (int ix=0; ix<nx && passed; ix++) {
      int i = ix + nx*(iy);
      if (a[i] != b[i]) passed = false;
    }
  }

  unit_func("read,write");
  unit_assert(passed);

  unit_finalize();
  PARALLEL_EXIT;

}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_FileHdf5.def.h)
