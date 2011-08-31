// See LICENSE_CELLO file for license and copyright information

/// @file     test_FileHdf5.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:47:35 PST 2008
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

  int nx = 17;
  int ny = 25;

  char * a_char = new char[nx*ny];
  char * b_char = new char[nx*ny];
  int * a_int = new int[nx*ny];
  int * b_int = new int[nx*ny];
  long * a_long = new long[nx*ny];
  long * b_long = new long[nx*ny];
  float * a_float = new float[nx*ny];
  float * b_float = new float[nx*ny];
  double * a_double = new double[nx*ny];
  double * b_double = new double[nx*ny];

  // Initialize a, clear b

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {
      int i = ix + nx*(iy);
      a_char[i]   = ix*3 + iy*5;
      a_int[i]    = ix*3 + iy*5;
      a_long[i]   = ix*3 + iy*5;
      a_float[i]  = ix*3 + iy*5;
      a_double[i] = ix*3 + iy*5;

      b_char[i]   = 0;
      b_int[i]    = 0;
      b_long[i]   = 0;
      b_float[i]  = 0;
      b_double[i] = 0;

    }
  }

  //--------------------------------------------------
  // Open file A for writing
  //--------------------------------------------------

   unit_func("open");

   int a_nx = nx;
   int a_ny = ny;
   int a_nz = 1;

   FileHdf5 hdf5_a("./","test_disk.h5","w");

   hdf5_a.open();

   hdf5_a.data_set   ("char",scalar_type_char, a_nx,a_ny,a_nz);
   hdf5_a.data_write (a_char);

   hdf5_a.data_set   ("int",scalar_type_int, a_nx,a_ny,a_nz);
   hdf5_a.data_write (a_int);

   hdf5_a.data_set   ("long",scalar_type_long, a_nx,a_ny,a_nz);
   hdf5_a.data_write (a_long);

   hdf5_a.data_set   ("float",scalar_type_float, a_nx,a_ny,a_nz);
   hdf5_a.data_write (a_float);

   hdf5_a.data_set   ("double",scalar_type_double, a_nx,a_ny,a_nz);
   hdf5_a.data_write (a_double);

   hdf5_a.close();

  //--------------------------------------------------
  // Open a file for reading
  //--------------------------------------------------

  int b_nx;
  int b_ny;
  int b_nz;
  scalar_type scalar_type;

  unit_func("data_get");

  FileHdf5 hdf5_b("./","test_disk.h5","r");
  hdf5_b.open();

  //----------------------------------------------------------------------
  unit_func("read double");
  //----------------------------------------------------------------------

  hdf5_b.data_get  ("double",&scalar_type, &b_nx,&b_ny,&b_nz);
  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);
  unit_assert (scalar_type == scalar_type_double);
  hdf5_b.data_read (b_double);

  bool p_double = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*(iy);
      if (a_double[i] != b_double[i]) p_double = false;
    }
  }

  unit_assert(p_double);

  //----------------------------------------------------------------------
  unit_func("read float");
  //----------------------------------------------------------------------

  hdf5_b.data_get   ("float",&scalar_type, &b_nx,&b_ny,&b_nz);
  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);
  unit_assert (scalar_type == scalar_type_float);
  hdf5_b.data_read (b_float);

  bool p_float  = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*(iy);
      if (a_float[i]  != b_float[i])  p_float  = false;
    }
  }

  unit_assert(p_float);

  //----------------------------------------------------------------------
  unit_func("read long");
  //----------------------------------------------------------------------

  hdf5_b.data_get   ("long",&scalar_type, &b_nx,&b_ny,&b_nz);
  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);
  unit_assert (scalar_type == scalar_type_long);
  hdf5_b.data_read (b_long);

  bool p_long   = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*(iy);
      if (a_long[i]   != b_long[i])   p_long   = false;
    }
  }

  unit_assert(p_long);

  //----------------------------------------------------------------------
  unit_func("read int");
  //----------------------------------------------------------------------

  hdf5_b.data_get   ("int",&scalar_type, &b_nx,&b_ny,&b_nz);
  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);
  unit_assert (scalar_type == scalar_type_int);
  hdf5_b.data_read (b_int);

  bool p_int    = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*(iy);
      if (a_int[i]    != b_int[i])    p_int    = false;
    }
  }

  unit_assert(p_int);

  //----------------------------------------------------------------------
  unit_func("read char");
  //----------------------------------------------------------------------

  hdf5_b.data_get   ("char",&scalar_type, &b_nx,&b_ny,&b_nz);
  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);
  unit_assert (scalar_type == scalar_type_char);
  hdf5_b.data_read (b_char);

  bool p_char   = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*(iy);
      if (a_char[i]   != b_char[i])   p_char   = false;
    }
  }

  unit_assert(p_char);


  hdf5_b.close();

  //--------------------------------------------------
  // Finalize
  //--------------------------------------------------

  delete [] a_char;
  delete [] a_int;
  delete [] a_long;
  delete [] a_float;
  delete [] a_double;
  delete [] b_char;
  delete [] b_int;
  delete [] b_long;
  delete [] b_float;
  delete [] b_double;

  unit_finalize();

  PARALLEL_EXIT;

}
PARALLEL_MAIN_END
