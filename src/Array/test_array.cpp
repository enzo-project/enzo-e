//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//
/** 
 *********************************************************************
 *
 * @file      test_array.cpp
 * @brief     Program implementing unit tests for the Array class
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      Thu Feb 21 16:04:03 PST 2008
 * @ingroup   Array
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include "cello.h"

#include "error.hpp"
#include "test.hpp"
#include "array.hpp"

int main()
{

  unit_class ("Array");

  unit_open();

  //----------------------------------------------------------------------
  // test single array WITH resize: length, size, clear, values, access
  //----------------------------------------------------------------------

  {
    Array a;

    int n0=10,n1=15,n2=20;
    a.resize(n0,n1,n2);
    int n = n0*n1*n2;
    int m0,m1,m2;
    a.size(&m0,&m1,&m2);
    int m = a.length();
    unit_func("length");
    unit_assert(n == m);
    unit_func("size");
    unit_assert(m0==n0);
    unit_assert(m1==n1);
    unit_assert(m2==n2);

    Scalar * av = a.values();

    // set values
    int i0,i1,i2;
    for (i0=0; i0<n0; i0++) {
      for (i1=0; i1<n1; i1++) {
	for (i2=0; i2<n2; i2++) {
	  int i = i0 + n0*(i1 + n1*i2);
	  av [ i ] = 100+i;
	}
      }
    }    

    // test values

    bool passed = true;
    for (i0=0; i0<n0; i0++) {
      for (i1=0; i1<n1; i1++) {
	for (i2=0; i2<n2; i2++) {
	  int i = i0 + n0*(i1 + n1*i2);
	  if (av[i] != a(i0,i1,i2) && passed) {
	    passed = false;
	    printf ("av[%d] = %g  a(%d,%d,%d) = %g\n",
		    i,av[i],i0,i1,i2,a(i0,i1,i2));
	  }
	}
      }
    }    
    unit_func("operator()");
    unit_assert(passed);

    // clear() to 0

    a.clear();

    passed = true;
    for (i0=0; i0<n0; i0++) {
      for (i1=0; i1<n1; i1++) {
	for (i2=0; i2<n2; i2++) {
	  int i = i0 + n0*(i1 + n1*i2);
	  if (av[i]!=0.0 && passed) {
	    passed = false;
	    printf ("av[%d] = %g\n",i,av[i]);
	  }
	}
      }
    }    
    unit_func("clear");
    unit_assert(passed);

  }
  //----------------------------------------------------------------------
  // test single array WITHOUT resize: length, size, clear, values, access
  //----------------------------------------------------------------------

  {

    int n0=10,n1=15,n2=20;
    
    Array a;
    a.resize(n0,n1,n2);
 
    int n = n0*n1*n2;
    int m0,m1,m2;
    a.size(&m0,&m1,&m2);
    int m = a.length();
    unit_func("length");
    unit_assert(n == m);
    unit_func("size");
    unit_assert(m0==n0);
    unit_assert(m1==n1);
    unit_assert(m2==n2);

    Scalar * av = a.values();

    // set values
    int i0,i1,i2;
    for (i0=0; i0<n0; i0++) {
      for (i1=0; i1<n1; i1++) {
	for (i2=0; i2<n2; i2++) {
	  int i = i0 + n0*(i1 + n1*i2);
	  av [ i ] = 100+i;
	}
      }
    }    

    // test values

    bool passed = true;
    for (i0=0; i0<n0; i0++) {
      for (i1=0; i1<n1; i1++) {
	for (i2=0; i2<n2; i2++) {
	  int i = i0 + n0*(i1 + n1*i2);
	  if (av[i]!=a(i0,i1,i2) && passed) {
	    passed = false;
	    printf ("av[%d] = %g  a(%d,%d,%d) = %g\n",
		    i,av[i],i0,i1,i2,a(i0,i1,i2));
	  }
	}
      }
    }    
    unit_func("operator()");
    unit_assert(passed);

    // Multiple arrays: copy
    Array b;
    b = a;

    Scalar * bv = b.values();
    
    // compare values

    passed = true;
    for (i0=0; i0<n0; i0++) {
      for (i1=0; i1<n1; i1++) {
	for (i2=0; i2<n2; i2++) {
	  int i = i0 + n0*(i1 + n1*i2);
	  if (bv[i]!=av[i] && passed) {
	    passed = false;
	    printf ("bv[%d] = %g  av[%d] = %g\n",i,bv[i],i,av[i]);
	  }
	}
      }
    }    
    unit_func("operator =");
    unit_assert(passed);

    // clear() to non-0

    a.clear(1.5);
    passed = true;
    for (i0=0; i0<n0; i0++) {
      for (i1=0; i1<n1; i1++) {
	for (i2=0; i2<n2; i2++) {
	  int i = i0 + n0*(i1 + n1*i2);
	  if (av[i]!=1.5 && passed) {
	    passed = false;
	    printf ("av[%d] = %g\n",i,av[i]);
	  }
	}
      }
    }    
    unit_func("clear");
    unit_assert(passed);

  }

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
  
  unit_close();

}
