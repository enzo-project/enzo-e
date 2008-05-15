/** 
*********************************************************************
*
* @file      time_matvec.cpp
* @brief     Program timing sparse matrix-vector multiplies
* @author    James Bordner
* @date      Thu Apr 24 16:25:43 PDT 2008
*
*********************************************************************
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "test_unit_.hpp"
#include "data_scalar_.hpp"
#include "data_array_.hpp"
#include "timer.hpp"

main(int argc, char ** argv)
{
  Timer timer_const;
  Timer overhead;
  Timer timer_symm;
  Timer timer_full;
  int i,j,k,i0,i1,i2;

  const int CACHE_SIZE = 4096;
  char cache[CACHE_SIZE];

  int kbegin,kend;

  if (argc >= 1) {
    kbegin = 8;
    kend   = 192;
  }
  if (argc >= 2) {
    kbegin = atoi(argv[1]);
    kend   = atoi(argv[1]);
  }
  if (argc >= 3) {
    kend   = atoi(argv[2]);
  }

  for (k=kbegin; k<=kend; k++) {

    int k2  = k + 2;
    int k22 = k2*k2;

    int count = 1000000/k/k/k+1;
    {
      Array A0(k2,k2,k2);  Scalar  * a0 = A0.values();
      Array AXP(k2,k2,k2); Scalar * axp = AXP.values();
      Array AXM(k2,k2,k2); Scalar * axm = AXM.values();
      Array AYP(k2,k2,k2); Scalar * ayp = AYP.values();
      Array AYM(k2,k2,k2); Scalar * aym = AYM.values();
      Array AZP(k2,k2,k2); Scalar * azp = AZP.values();
      Array AZM(k2,k2,k2); Scalar * azm = AZM.values();
      Array   X(k2,k2,k2); Scalar   * x =   X.values();
      Array   B(k2,k2,k2); Scalar    *b =   B.values();

      if (a0 && axp && axm && ayp && aym && azp && azm) {

	for (j=0; j<count; j++) {

	  timer_const.start();

	  for (i0=1; i0<k+1; i0++) {
	    for (i1=1; i1<k+1; i1++) {
	      for (i2=1; i2<k+1; i2++) {
		i = i0 + k2*(i1 + k2*i2);
		b[i] = azm[0]*x[i-k22]
		  +    aym[0]*x[i-k2]
		  +    axm[0]*x[i-1]
		  +     a0[0]*x[i]
		  +    axp[0]*x[i+1]
		  +    ayp[0]*x[i+k2]
		  +    azp[0]*x[i+k22];
	      }
	    }
	  }

	  timer_const.stop();

	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}

	for (j=0; j<count; j++) {
	  overhead.start();
	  overhead.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}
      }
    }

    {
      Array A0(k2,k2,k2);  Scalar  * a0 = A0.values();
      Array AXP(k2,k2,k2); Scalar * axp = AXP.values();
      Array AXM(k2,k2,k2); Scalar * axm = AXM.values();
      Array AYP(k2,k2,k2); Scalar * ayp = AYP.values();
      Array AYM(k2,k2,k2); Scalar * aym = AYM.values();
      Array AZP(k2,k2,k2); Scalar * azp = AZP.values();
      Array AZM(k2,k2,k2); Scalar * azm = AZM.values();
      Array   X(k2,k2,k2); Scalar   * x =   X.values();
      Array   B(k2,k2,k2); Scalar    *b =   B.values();
    
      if (a0 && axp && axm && ayp && aym && azp && azm) {


	for (j=0; j<count; j++) {
	  timer_symm.start();
	  for (i0=1; i0<k+1; i0++) {
	    for (i1=1; i1<k+1; i1++) {
	      for (i2=1; i2<k+1; i2++) {
		i = i0 + k2*(i1 + k2*i2);
		b[i] = azp[i]*x[i-k22]
		  +    ayp[i]   *x[i-k2]
		  +    axp[i]    *x[i-1]
		  +     a0[i]      *x[i]
		  +    axp[i+1]    *x[i+1]
		  +    ayp[i+k2]   *x[i+k2]
		  +    azp[i+k22]*x[i+k22];
	      }
	    }
	  }
	  timer_symm.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}
      }
    }

    {
      Array A0(k2,k2,k2);  Scalar  * a0 = A0.values();
      Array AXP(k2,k2,k2); Scalar * axp = AXP.values();
      Array AXM(k2,k2,k2); Scalar * axm = AXM.values();
      Array AYP(k2,k2,k2); Scalar * ayp = AYP.values();
      Array AYM(k2,k2,k2); Scalar * aym = AYM.values();
      Array AZP(k2,k2,k2); Scalar * azp = AZP.values();
      Array AZM(k2,k2,k2); Scalar * azm = AZM.values();
      Array   X(k2,k2,k2); Scalar   * x =   X.values();
      Array   B(k2,k2,k2); Scalar    *b =   B.values();
    
      if (a0 && axp && axm && ayp && aym && azp && azm) {

	for (j=0; j<count; j++) {
	timer_full.start();
	  for (i0=1; i0<k+1; i0++) {
	    for (i1=1; i1<k+1; i1++) {
	      for (i2=1; i2<k+1; i2++) {
		i = i0 + k2*(i1 + k2*i2);
		b[i] = azm[i-k22]*x[i-k22]
		  +    aym[i-k2]   *x[i-k2]
		  +    axm[i-1]    *x[i-1]
		  +     a0[i]      *x[i]
		  +    axp[i+1]    *x[i+1]
		  +    ayp[i+k2]   *x[i+k2]
		  +    azp[i+k22]*x[i+k22];
	      }
	    }
	  }
	  timer_full.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}
      }
    }


    int num_flops = count*k*k*k*13;
    printf ("%d %20.12f %20.12f %20.12f \n",
	    k,
	    1e-6*num_flops/(timer_const.value()-overhead.value()),
	    1e-6*num_flops/(timer_symm.value()-overhead.value()),
	    1e-6*num_flops/(timer_full.value()-overhead.value())
	    ); 
    fflush(stdout);

    timer_const.clear();
    timer_symm.clear();
    timer_full.clear();
    overhead.clear();

  } // for
}
