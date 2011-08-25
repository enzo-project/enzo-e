// See LICENSE_CELLO file for license and copyright information

/// @file     test_FileHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:47:35 PST 2008
/// @todo     Multiple datasets with different precisions
/// @brief    Program implementing unit tests for the FileHdf5 class
 
#include "main.hpp" 
#include "test.hpp"

#include "disk.hpp"

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;

  unit_init();

  unit_class("FileHdf5");

  //--------------------------------------------------
  // Initialize
  //--------------------------------------------------

  // Allocate arrays
  int nx = 100;
  int ny = 50;
  double * a = new double[nx*ny];
  double * b = new double[nx*ny];

  // Initialize array "a" only

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {
      int i = ix + nx*(iy);
      a[i] = ix*3 + iy*5;
    }
  }

  FileHdf5 hdf5;

  //--------------------------------------------------
  // Open a file for writing
  //--------------------------------------------------

  unit_func("open");

  int mx,my,mz;
  hdf5.open("open_test.h5","w");

  unit_assert(true);

  // Open a dataset
  mx = nx;
  my = ny;
  mz = 1;
  hdf5.open_dataset ("dataset",precision_double,mx,my,mz);

  // Write the dataset
  hdf5.write((char *)a,precision_double);

  // Close the dataset
  hdf5.close_dataset ();

  // Close the file
  hdf5.close();

  //--------------------------------------------------
  // Open a file for reading
  //--------------------------------------------------

  hdf5.open("open_test.h5","r");

  // Open a dataset
  hdf5.open_dataset ("dataset",&mx,&my,&mz);

  // Read the dataset
  hdf5.read((char *)b,precision_double);

  // Close the dataset
  hdf5.close_dataset ();

  // Close the file
  hdf5.close();

  //--------------------------------------------------
  // Compare data written to data read
  //--------------------------------------------------

  bool passed = true;

  // Compare a with b
  for (int iy=0; iy<ny && passed; iy++) {
    for (int ix=0; ix<nx && passed; ix++) {
      int i = ix + nx*(iy);
      if (a[i] != b[i]) passed = false;
    }
  }

  // Report test results
  unit_func("read");
  unit_assert(passed);
  unit_func("write");
  unit_assert(passed);

  //--------------------------------------------------
  // Finalize
  //--------------------------------------------------

  delete [] a;
  delete [] b;

  unit_finalize();

  PARALLEL_EXIT;

}
PARALLEL_MAIN_END
