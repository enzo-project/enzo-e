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

  unit_init(0,1);

  unit_class("FileHdf5");

  //--------------------------------------------------
  // Initialize
  //--------------------------------------------------

  // data arrays

  const int nx = 7;
  const int ny = 5;

  char      * a_char      = new char      [nx*ny];
  int       * a_int       = new int       [nx*ny];
  long      * a_long      = new long      [nx*ny];
  long long * a_long_long = new long long [nx*ny];
  float     * a_float     = new float     [nx*ny];
  double    * a_double    = new double    [nx*ny];

  char      * b_char      = new char      [nx*ny];
  int       * b_int       = new int       [nx*ny];
  long      * b_long      = new long      [nx*ny];
  long long * b_long_long = new long long [nx*ny];
  float     * b_float     = new float     [nx*ny];
  double    * b_double    = new double    [nx*ny];

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {

      int i = ix + nx*iy;

      a_char[i]   = ix*3 + iy*5;
      a_int[i]    = ix*3 + iy*5;
      a_long[i]   = ix*3 + iy*5;
      a_long_long[i]   = ix*3 + iy*5;
      a_float[i]  = ix*3 + iy*5;
      a_double[i] = ix*3 + iy*5;

      b_char[i]   = 0;
      b_int[i]    = 0;
      b_long[i]   = 0;
      b_long_long[i]   = 0;
      b_float[i]  = 0;
      b_double[i] = 0;

    }
  }

  // Initialize meta arrays

  const int mx = 3;
  const int my = 2;

  char ma_char   [mx*my];
  int  ma_int    [mx*my];
  long ma_long   [mx*my];
  long ma_long_long   [mx*my];
  float ma_float  [mx*my];
  double ma_double[mx*my];

  char mb_char   [mx*my];
  int mb_int    [mx*my];
  long mb_long   [mx*my];
  long mb_long_long   [mx*my];
  float mb_float  [mx*my];
  double mb_double[mx*my];

  for (int iy=0; iy<my; iy++) {
    for (int ix=0; ix<mx; ix++) {
      int i = ix + mx*iy;
      ma_char[i]   = ix*5 - iy*13;
      ma_int[i]    = ix*5 - iy*13;
      ma_long[i]   = ix*5 - iy*13;
      ma_long_long[i]   = ix*5 - iy*13;
      ma_float[i]  = ix*5 - iy*13;
      ma_double[i] = ix*5 - iy*13;

      mb_char[i]   = 0;
      mb_int[i]    = 0;
      mb_long[i]   = 0;
      mb_long_long[i]   = 0;
      mb_float[i]  = 0;
      mb_double[i] = 0;

    }
  }

  //--------------------------------------------------
  // Create a file
  //--------------------------------------------------

  unit_func("open()");

  int a_nx = nx;
  int a_ny = ny;
  int a_nz = 1;

  int a_mx = mx;
  int a_my = my;

  FileHdf5 hdf5_a("./","test_disk.h5");

  hdf5_a.file_create();

  hdf5_a.file_write_meta(&mx, "mx", scalar_type_int);
  hdf5_a.file_write_meta(&my, "my", scalar_type_int);

  hdf5_a.group_create ("/test");
  hdf5_a.data_create ("char", scalar_type_char, a_nx,a_ny,a_nz);
  hdf5_a.data_write_meta(&ma_char, "meta_char", scalar_type_char, a_mx, a_my);
  hdf5_a.data_write (a_char);
  hdf5_a.data_close ();
  hdf5_a.group_close ();

  hdf5_a.group_create ("/test/int");
  hdf5_a.data_create ("int", scalar_type_int, a_nx,a_ny,a_nz);
  hdf5_a.data_write_meta(&ma_int, "meta_int", scalar_type_int, a_mx, a_my);
  hdf5_a.data_write (a_int);
  hdf5_a.data_close ();
  hdf5_a.group_close();

  hdf5_a.group_create ("/test2/long/long/type/is/long");

  hdf5_a.data_create ("long_long", scalar_type_long_long, a_nx,a_ny,a_nz);
  hdf5_a.data_write_meta
    (&ma_long_long, "meta_long_long", scalar_type_long_long, a_mx, a_my);
  hdf5_a.data_write (a_long_long);
  hdf5_a.data_close ();
  hdf5_a.group_close();

  hdf5_a.group_create ("/test2/scalar/long/group");

  hdf5_a.data_create ("long", scalar_type_long, a_nx,a_ny,a_nz);
  hdf5_a.data_write_meta(&ma_long, "meta_long", scalar_type_long, a_mx, a_my);
  hdf5_a.data_write (a_long);
  hdf5_a.data_close ();

  hdf5_a.data_create ("float", scalar_type_float, a_nx,a_ny,a_nz);
  hdf5_a.data_write_meta(&ma_float, "meta_float", scalar_type_float, a_mx, a_my);
  hdf5_a.data_write (a_float);
  hdf5_a.data_close ();
  hdf5_a.group_close();

  hdf5_a.data_create ("double",scalar_type_double, a_nx,a_ny,a_nz);
  hdf5_a.data_write_meta(&ma_double, "meta_double", scalar_type_double, a_mx, a_my);
  hdf5_a.data_write (a_double);
  hdf5_a.data_close ();


  hdf5_a.file_close();

  //--------------------------------------------------
  // Reopen the file
  //--------------------------------------------------

  int b_nx;
  int b_ny;
  int b_nz;
  int b_mx;
  int b_my;

  scalar_type type;

  unit_func("file_open()");

  FileHdf5 hdf5_b("./","test_disk.h5");
  hdf5_b.file_open();

  unit_assert(true);
 
  //----------------------------------------------------------------------
  unit_func("file_read_meta()");
  //----------------------------------------------------------------------

  type = scalar_type_unknown;
  hdf5_b.file_read_meta(&b_mx, "mx", &type);

  unit_assert (a_mx == b_mx);
  unit_assert (type == scalar_type_int);
  type = scalar_type_unknown;
  hdf5_b.file_read_meta(&b_my, "my", &type);

  unit_assert (a_my == b_my);
  unit_assert (type == scalar_type_int);

  //======================================================================
  unit_func("double data_open()");
  //======================================================================

  type = scalar_type_unknown;
  hdf5_b.data_open ("double",&type, &b_nx,&b_ny,&b_nz);

  unit_assert (type == scalar_type_double);

  //----------------------------------------------------------------------
  unit_func("double data_read_meta()");
  //----------------------------------------------------------------------

  type = scalar_type_unknown;
  hdf5_b.data_read_meta (&mb_double, "meta_double", &type, &b_mx, &b_my);

  unit_assert (a_mx == b_mx);
  unit_assert (a_my == b_my);
  unit_assert (type == scalar_type_double);

  //----------------------------------------------------------------------
  unit_func("double data_read()");
  //----------------------------------------------------------------------

  hdf5_b.data_read (b_double);

  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);

  //----------------------------------------------------------------------
  unit_func("double data match");
  //----------------------------------------------------------------------

  bool p_double = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*iy;
      if ( ! (a_double[i] == b_double[i]) ) {
	printf ("MISMATCH double meta %d %d  %g %g\n",ix,iy,a_double[i],b_double[i]);
      }
      p_double = p_double && (a_double[i] == b_double[i]);
    }
  }

  unit_assert(p_double);

  //----------------------------------------------------------------------
  unit_func("double meta match");
  //----------------------------------------------------------------------

  bool mp_double = true;

  for (int iy=0; iy<my ; iy++) {
    for (int ix=0; ix<mx ; ix++) {
      int i = ix + mx*iy;
      mp_double = mp_double && (ma_double[i] == mb_double[i]);
    }
  }

  unit_assert(mp_double);

  hdf5_b.data_close ();

  //======================================================================
  unit_func("char data_open()");
  //======================================================================

  hdf5_b.group_open ("/test");

  type = scalar_type_unknown;
  hdf5_b.data_open ("char",&type, &b_nx,&b_ny,&b_nz);

  unit_assert (type == scalar_type_char);

  //----------------------------------------------------------------------
  unit_func("char data_read_meta()");
  //----------------------------------------------------------------------

  type = scalar_type_unknown;
  hdf5_b.data_read_meta (&mb_char, "meta_char", &type, &b_mx, &b_my);

  unit_assert (a_mx == b_mx);
  unit_assert (a_my == b_my);
  unit_assert (type == scalar_type_char);

  //----------------------------------------------------------------------
  unit_func("char data_read()");
  //----------------------------------------------------------------------

  hdf5_b.data_read (b_char);

  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);

  //----------------------------------------------------------------------
  unit_func("char data match");
  //----------------------------------------------------------------------

  bool p_char = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*iy;
      p_char = p_char && (a_char[i] == b_char[i]);
    }
  }

  unit_assert(p_char);

  //----------------------------------------------------------------------
  unit_func("char meta match");
  //----------------------------------------------------------------------

  bool mp_char = true;

  for (int iy=0; iy<my ; iy++) {
    for (int ix=0; ix<mx ; ix++) {
      int i = ix + mx*iy;
      if ( ! (ma_char[i] == mb_char[i]) ) {
	printf ("MISMATCH char meta %d %d  %d %d\n",ix,iy,ma_char[i],mb_char[i]);
      }
      mp_char = mp_char && (ma_char[i] == mb_char[i]);
    }
  }

  unit_assert(mp_char);

  hdf5_b.data_close ();

  hdf5_b.group_close();

  //======================================================================
  unit_func("int data_open()");
  //======================================================================

  hdf5_b.group_open ("/test/int");

  type = scalar_type_unknown;
  hdf5_b.data_open ("int",&type, &b_nx,&b_ny,&b_nz);

  unit_assert (type == scalar_type_int);

  //----------------------------------------------------------------------
  unit_func("int data_read_meta()");
  //----------------------------------------------------------------------

  type = scalar_type_unknown;
  hdf5_b.data_read_meta (&mb_int, "meta_int", &type, &b_mx, &b_my);

  unit_assert (a_mx == b_mx);
  unit_assert (a_my == b_my);
  unit_assert (type == scalar_type_int);

  //----------------------------------------------------------------------
  unit_func("int data_read()");
  //----------------------------------------------------------------------

  hdf5_b.data_read (b_int);

  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);

  //----------------------------------------------------------------------
  unit_func("int data match");
  //----------------------------------------------------------------------

  bool p_int = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*iy;
      p_int = p_int && (a_int[i] == b_int[i]);
    }
  }

  unit_assert(p_int);

  //----------------------------------------------------------------------
  unit_func("int meta match");
  //----------------------------------------------------------------------

  bool mp_int = true;

  for (int iy=0; iy<my ; iy++) {
    for (int ix=0; ix<mx ; ix++) {
      int i = ix + mx*iy;
      if ( ! (ma_int[i] == mb_int[i]) ) {
	printf ("MISMATCH int meta %d %d  %d %d\n",ix,iy,ma_int[i],mb_int[i]);
      }
      mp_int = mp_int && (ma_int[i] == mb_int[i]);
    }
  }

  unit_assert(mp_int);

  hdf5_b.data_close ();
  hdf5_b.group_close ();

  //======================================================================
  unit_func("float data_open()");
  //======================================================================

  hdf5_b.group_open ("/test2/scalar/long/group");

  type = scalar_type_unknown;
  hdf5_b.data_open ("float",&type, &b_nx,&b_ny,&b_nz);

  unit_assert (type == scalar_type_float);

  //----------------------------------------------------------------------
  unit_func("float data_read_meta()");
  //----------------------------------------------------------------------

  type = scalar_type_unknown;
  hdf5_b.data_read_meta (&mb_float, "meta_float", &type, &b_mx, &b_my);

  unit_assert (a_mx == b_mx);
  unit_assert (a_my == b_my);
  unit_assert (type == scalar_type_float);

  //----------------------------------------------------------------------
  unit_func("float data_read()");
  //----------------------------------------------------------------------

  hdf5_b.data_read (b_float);

  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);

  //----------------------------------------------------------------------
  unit_func("float data match");
  //----------------------------------------------------------------------

  bool p_float = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*iy;
      p_float = p_float && (a_float[i] == b_float[i]);
    }
  }

  unit_assert(p_float);

  hdf5_b.data_close ();

  hdf5_b.group_close();

  //----------------------------------------------------------------------
  unit_func("float meta match");
  //----------------------------------------------------------------------

  bool mp_float = true;

  for (int iy=0; iy<my ; iy++) {
    for (int ix=0; ix<mx ; ix++) {
      int i = ix + mx*iy;
      if ( ! (ma_float[i] == mb_float[i]) ) {
	printf ("MISMATCH float meta %d %d  %g %g\n",ix,iy,ma_float[i],mb_float[i]);
      }
      mp_float = mp_float && (ma_float[i] == mb_float[i]);
    }
  }

  unit_assert(mp_float);

  //======================================================================
  unit_func("long data_open()");
  //======================================================================

  hdf5_b.group_open ("/test2/scalar/long/group");
  
  type = scalar_type_unknown;
  hdf5_b.data_open ("long",&type, &b_nx,&b_ny,&b_nz);

  unit_assert (type == 
	       ((sizeof(int)==sizeof(long)) ? 
		scalar_type_int : scalar_type_long));

  //----------------------------------------------------------------------
  unit_func("long data_read_meta()");
  //----------------------------------------------------------------------

  type = scalar_type_unknown;
  hdf5_b.data_read_meta (&mb_long, "meta_long", &type, &b_mx, &b_my);

  unit_assert (a_mx == b_mx);
  unit_assert (a_my == b_my);
  unit_assert (type == 
	       ((sizeof(int)==sizeof(long)) ? 
		scalar_type_int : scalar_type_long));

  //----------------------------------------------------------------------
  unit_func("long data_read()");
  //----------------------------------------------------------------------

  hdf5_b.data_read (b_long);

  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);

  //----------------------------------------------------------------------
  unit_func("long data match");
  //----------------------------------------------------------------------

  bool p_long = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*iy;
      p_long = p_long && (a_long[i] == b_long[i]);
    }
  }

  unit_assert(p_long);

  hdf5_b.data_close ();

  hdf5_b.group_close ();

  //----------------------------------------------------------------------
  unit_func("long meta match");
  //----------------------------------------------------------------------

  bool mp_long = true;

  for (int iy=0; iy<my ; iy++) {
    for (int ix=0; ix<mx ; ix++) {
      int i = ix + mx*iy;
      if ( ! (ma_long[i] == mb_long[i]) ) {
	printf ("MISMATCH long meta %d %d  %ld %ld\n",ix,iy,ma_long[i],mb_long[i]);
      }
      mp_long = mp_long && (ma_long[i] == mb_long[i]);
    }
  }

  unit_assert(mp_long);

  //======================================================================
  unit_func("long long data_open()");
  //======================================================================

  hdf5_b.group_open ("/test2/long/long/type/is/long");
  
  type = scalar_type_unknown;
  hdf5_b.data_open ("long_long",&type, &b_nx,&b_ny,&b_nz);
  scalar_type scalar_type_expected = 
    ((sizeof(long)==sizeof(long long)) ? 
     scalar_type_long : scalar_type_long_long);

  unit_assert (type == scalar_type_expected);

  //----------------------------------------------------------------------
  unit_func("long long data_read_meta()");
  //----------------------------------------------------------------------

  type = scalar_type_unknown;
  hdf5_b.data_read_meta (&mb_long_long, "meta_long_long", &type, &b_mx, &b_my);

  unit_assert (a_mx == b_mx);
  unit_assert (a_my == b_my);
  unit_assert (type == ((sizeof(long)==sizeof(long long)) ? scalar_type_long : scalar_type_long_long));

  //----------------------------------------------------------------------
  unit_func("long long data_read()");
  //----------------------------------------------------------------------

  hdf5_b.data_read (b_long_long);

  unit_assert (a_nx == b_nx);
  unit_assert (a_ny == b_ny);
  unit_assert (a_nz == b_nz);

  //----------------------------------------------------------------------
  unit_func("long long data match");
  //----------------------------------------------------------------------

  bool p_long_long = true;

  for (int iy=0; iy<ny ; iy++) {
    for (int ix=0; ix<nx ; ix++) {
      int i = ix + nx*iy;
      p_long_long = p_long_long && (a_long_long[i] == b_long_long[i]);
    }
  }

  unit_assert(p_long_long);

  hdf5_b.data_close ();
  hdf5_b.group_close ();

  //----------------------------------------------------------------------
  unit_func("long long meta match");
  //----------------------------------------------------------------------

  bool mp_long_long = true;

  for (int iy=0; iy<my ; iy++) {
    for (int ix=0; ix<mx ; ix++) {
      int i = ix + mx*iy;
      if ( ! (ma_long_long[i] == mb_long_long[i]) ) {
	printf ("MISMATCH long long meta %d %d  %ld %ld\n",ix,iy,ma_long_long[i],mb_long_long[i]);
      }
      mp_long_long = mp_long_long && (ma_long_long[i] == mb_long_long[i]);
    }
  }

  unit_assert(mp_long_long);

  hdf5_b.file_close();

  //--------------------------------------------------
  // Finalize
  //--------------------------------------------------

  delete [] a_char;
  delete [] a_int;
  delete [] a_long;
  delete [] a_long_long;
  delete [] a_float;
  delete [] a_double;
  delete [] b_char;
  delete [] b_int;
  delete [] b_long;
  delete [] b_long_long;
  delete [] b_float;
  delete [] b_double;

  unit_finalize();

 PARALLEL_EXIT;

}
PARALLEL_MAIN_END
