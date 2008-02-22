/** 
 *********************************************************************
 *
 * @file      test_array.cpp
 * @brief     Program implementing unit tests for the Array class
 * @author    James Bordner
 * @date      Thu Feb 21 16:04:03 PST 2008
 *
 *********************************************************************
 */
 
#include <stdio.h>

#include "scalar.hpp"
#include "array.hpp"
#include "unit.hpp"

main()
{
  UNIT_CLASS ("Array");

  //----------------------------------------------------------------------
  // test single array with resize: length, size, and values, and element access
  //----------------------------------------------------------------------

  {
    Array a;
 
    int n0=10,n1=15,n2=20;
    a.resize(n0,n1,n2);
    int n = n0*n1*n2;
    int m0,m1,m2,m3;
    a.size(&m0,&m1,&m2,&m3);
    int m = a.length();
    UNIT_FUNC("length");
    UNIT_ASSERT(n == m);
    UNIT_FUNC("size");
    UNIT_ASSERT(m0==n0);
    UNIT_ASSERT(m1==n1);
    UNIT_ASSERT(m2==n2);
    UNIT_ASSERT(m3==1);

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
    UNIT_FUNC("operator()");
    UNIT_ASSERT(passed);

  }
  //----------------------------------------------------------------------
  // test single array with resize: length, size, and values, and element access
  //----------------------------------------------------------------------

  {

    int n0=10,n1=15,n2=20;
    
    Array a(n0,n1,n2);
 
    int n = n0*n1*n2;
    int m0,m1,m2,m3;
    a.size(&m0,&m1,&m2,&m3);
    int m = a.length();
    UNIT_FUNC("length");
    UNIT_ASSERT(n == m);
    UNIT_FUNC("size");
    UNIT_ASSERT(m0==n0);
    UNIT_ASSERT(m1==n1);
    UNIT_ASSERT(m2==n2);
    UNIT_ASSERT(m3==1);

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
    UNIT_FUNC("operator()");
    UNIT_ASSERT(passed);

    // Multiple arrays: copy
    Array b;
    b.copy(a);

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
    UNIT_FUNC("copy");
    UNIT_ASSERT(passed);
  }

}
