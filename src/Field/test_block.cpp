// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the Block class

#include "cello.h"

#include "error.hpp"
#include "test.hpp"
#include "field.hpp"

int main()
{

  unit_class ("Block");

  {

    int n0=10,n1=15,n2=20;
    Array a;
    int m0,m1,m2;
    int nx,ny,nz,na;
    int mx,my,mz,ma;
    int i,j,k;
    bool is_same = true;

    a.resize(n0,n1,n2);

    a.size(&m0,&m1,&m2);
    Block b (a.values(),NULL,m0,m1,m2);

    unit_func("get_size");
    b.get_size(&nx,&ny,&nz,&na);
    printf ("(nx ny nz na) = %d %d %d %d\n",nx,ny,nz,na);
    unit_assert (nx == 10 && ny == 15 && nz == 20);

    Scalar * bv = b.values(0);
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  i = ix + nx*(iy + ny*iz);
	  bv[i] = ix + 7*iy + 19*iz;
	}
      }
    }

    unit_func("get_inc");
    b.get_inc(&mx,&my,&mz,&ma);
    printf ("(mx my mz ma) = %d %d %d %d\n",mx,my,mz,ma);
    unit_assert (mx == 1 && my == m0 && mz == m0*m1);
    is_same = true;
    i = 0;
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  j = ix + nx*(iy + ny*iz);
	  if (i != j) is_same = false;
	  k = ix*mx + iy*my + iz*mz;
	  if (j != k) is_same = false;
	  bv[i] = ix + 7*iy + 19*iz;
	  i += mx;
	}
	i += my - nx*mx;
      }
      i += mz - ny*my;
    }
    unit_assert (is_same);

    // Test treating as two interleaved arrays
    m0 /= 2;
    int p[] = {3,0,1,2};
    Block b2 (a.values(),p,m0,m1,m2,2);

    unit_func("get_size");
    b2.get_size(&nx,&ny,&nz,&na);
    printf ("(nx ny nz na) = %d %d %d %d\n",nx,ny,nz,na);
    unit_assert (nx == m0 && ny == m1 && nz == m2);

    bv = b2.values(0);
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  i = na*(ix + nx*(iy + ny*iz));
	  bv[i] = ix + 7*iy + 19*iz;
	}
      }
    }

    unit_func("get_inc");
    b2.get_inc(&mx,&my,&mz,&ma);
    printf ("(mx my mz ma) = %d %d %d %d\n",mx,my,mz,ma);
    unit_assert (mx == na && my == na*nx && mz == na*nx*ny);
    i = 0;
    is_same = true;
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int j = na*(ix + nx*(iy + ny*iz));
	  if (i != j) is_same = false;
	  k = ix*mx + iy*my + iz*mz;
	  if (i != k) is_same = false;
	  bv[i] = ix + 7*iy + 19*iz;
	  i += mx;
	}
	i += my - nx*mx;
      }
      i += mz-ny*my;
    }
    unit_assert (is_same);
  }
  
}
