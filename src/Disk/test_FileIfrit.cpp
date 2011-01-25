// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_FileIfrit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:47:35 PST 2008
/// @brief    Program implementing unit tests for the FileIfrit class
 
#include "cello.hpp"

#include "disk.hpp"

#include "test.hpp"

#include PARALLEL_CHARM_INCLUDE(test_FileIfrit.decl.h)

PARALLEL_MAIN_BEGIN

{
  PARALLEL_INIT;

  unit_init();

  int n0 = 64;
  int n1 = 64;
  int n2 = 64;

  const char filename[] = "FileIfrit_test.bin";
  int n = n0*n1*n2;

  float * a = new float[n];

  for (int i2=0; i2<n2; i2++) {
    double x2 = 1.0*i2/(n2-1);
    for (int i1=0; i1<n1; i1++) {
      double x1 = 1.0*i1/(n1-1);
      for (int i0=0; i0<n0; i0++) {
	double x0 = 1.0*i0/(n0-1);
	int i = i0 + n0*(i1 + n1*i2);
	a[i] = (x0-0.5)*(x0-0.5) + (x1-0.5)*(x1-0.5) + (x2-0.5)*(x2-0.5);
      }
    }
  }

  unit_class ("FileIfrit");

  FileIfrit ifrit;

  unit_func("write_bin");
  ifrit.write_bin(filename,a,n0,n1,n2);
  unit_assert(true);

  unit_func("read_bin");
  float * b = new float[n];
  int m0,m1,m2;
  ifrit.read_bin(filename,b,&m0,&m1,&m2);
  unit_assert (m0 == n0);
  unit_assert (m1 == n1);
  unit_assert (m2 == n2);

  bool passed = true;
  for (int i2=0; i2<n2 && passed; i2++) {
    for (int i1=0; i1<n1 && passed; i1++) {
      for (int i0=0; i0<n0 && passed; i0++) {
	int i = i0 + n0*(i1 + n1*i2);
	if (a[i] != b[i]) passed = false;
      }
    }
  }

  unit_assert(passed);

  unit_finalize();

  PARALLEL_EXIT;
}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_FileIfrit.def.h)
