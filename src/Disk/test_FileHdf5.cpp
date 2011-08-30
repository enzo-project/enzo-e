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
  double * a_double = new double[nx*ny];
  double * b_double = new double[nx*ny];

  // Initialize a, clear b

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {
      int i = ix + nx*(iy);
      a_double[i] = ix*3 + iy*5;
      b_double[i] = 0;
    }
  }

  //--------------------------------------------------
  // Open a file for writing
  //--------------------------------------------------

  unit_func("open");

  int a_nx, a_ny, a_nz;

  // Open a dataset
  a_nx = nx;
  a_ny = ny;
  a_nz = 1;

  {
    FileHdf5 hdf5("./","disk_double.h5","w");

    hdf5.open();

    // Create 1D array of integers
    hdf5.data_set ("double",scalar_double, a_nx,a_ny,a_nz);

    // Write the data
    hdf5.data_write (a_double);

    // Close the file
    hdf5.close();
  }

  //--------------------------------------------------
  // Open a file for reading
  //--------------------------------------------------

  {
    FileHdf5 hdf5("./","disk_double.h5","r");

    hdf5.open();

    // Get the dataset size and type

    unit_func("data_get");

    scalar_type b_type;
    int b_nx,b_ny,b_nz;
    hdf5.data_get ("double",&b_type, &b_nx, &b_ny, &b_nz);

    unit_assert (b_type == a_type);
    unit_assert (b_nx == a_nx);
    unit_assert (b_ny == a_ny);
    unit_assert (b_nz == a_nz);

    // Read the dataset
    hdf5.read((char *)b,precision_double);

    // Close the dataset
    hdf5.close_data ();

    // Close the file
    hdf5.close();

  }

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
