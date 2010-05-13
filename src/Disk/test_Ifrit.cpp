// $Id: test_Ifrit.cpp 1451 2010-05-10 18:23:30Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Ifrit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:47:35 PST 2008
/// @brief    Program implementing unit tests for the Ifrit class
 
#include <stdio.h>
#include <string>

#include "error.hpp"
#include "test.hpp"
#include "disk.hpp"

int main(int argc, char ** argv)
{
  int n0 = 64;
  int n1 = 64;
  int n2 = 64;

  int n = n0*n1*n2;

  float * a = new float[n];

  for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      for (int i0=0; i0<n0; i0++) {
	int i = i0 + n0*(i1 + n1*i2);
	a[i] = i0*3 + i1*5;
      }
    }
  }

  unit_class ("Ifrit");

  Ifrit ifrit;

  unit_func("write_bin");
  ifrit.write_bin("file_write_test.ifrit",a,n0,n1,n2);
  unit_assert(true);

  unit_func("read_bin");
  float * b = new float[n];
  int m0,m1,m2;
  ifrit.read_bin("file_write_test.ifrit",b,&m0,&m1,&m2);
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

}
