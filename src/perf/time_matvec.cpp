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

#define MIN(a,b) ((a)<(b)?(a):(b))
#define Scalar double

#include "timer.hpp"

main(int argc, char ** argv)
{
  Timer time_overhead;
  Timer time_const;
  Timer time_symm;
  Timer time_full;
  Timer time_block;
  int i,j,k,i0,i1,i2,ii2;

  const int CACHE_SIZE = 65536;
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

  Scalar * array = new Scalar [(kend+2)*(kend+2)*(kend+2)*9];

  for (k=kbegin; k<=kend; k++) {

    int k2  = k + 2;
    int k22 = k2*k2;
    int n   = k2*k2*k2;

    int count = 1000000/k/k/k + 1;

    int ind = 0;
    Scalar   * x = &array[ind]; ind += n;
    Scalar    *b = &array[ind]; ind += n;
    Scalar * azm = &array[ind]; ind += n;
    Scalar * aym = &array[ind]; ind += n;
    Scalar * axm = &array[ind]; ind += n;
    Scalar  * a0 = &array[ind]; ind += n;
    Scalar * axp = &array[ind]; ind += n;
    Scalar * ayp = &array[ind]; ind += n;
    Scalar * azp = &array[ind]; ind += n;

    {
      if (a0 && axp && axm && ayp && aym && azp && azm) {

	for (j=0; j<count; j++) {

	  time_const.start();

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

	  time_const.stop();

	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}

	for (j=0; j<count; j++) {
	  time_overhead.start();
	  time_overhead.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}
      }
    }

    {
    
      if (a0 && axp && axm && ayp && aym && azp && azm) {


	for (j=0; j<count; j++) {
	  time_symm.start();
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
	  time_symm.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}
      }
    }

    {
    
      if (a0 && axp && axm && ayp && aym && azp && azm) {

	for (j=0; j<count; j++) {
	  time_full.start();
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
	  time_full.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}
      }
    }

    {
      if (a0 && axp && axm && ayp && aym && azp && azm) {

	for (j=0; j<count; j++) {

	  time_block.start();

	  for (ii2=1; ii2<k+1; ii2+=4) {
	    for (i0=1; i0<k+1; i0++) {
	      for (i1=1; i1<k+1; i1++) {
		for (i2=ii2; i2<MIN(k+1,ii2+4); i2++) {
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
	  }

	  time_block.stop();

	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}

	for (j=0; j<count; j++) {
	  time_overhead.start();
	  time_overhead.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 1;
	  }
	}
      }
    }


    int num_flops = count*k*k*k*13;
    printf ("%d %20.12f %20.12f %20.12f %20.12f\n",
	    k,
	    1e-6*num_flops/(time_const.value()-time_overhead.value()),
	    1e-6*num_flops/(time_symm.value()-time_overhead.value()),
	    1e-6*num_flops/(time_full.value()-time_overhead.value()),
	    1e-6*num_flops/(time_block.value()-time_overhead.value())
	    ); 
    fflush(stdout);

    time_const.clear();
    time_symm.clear();
    time_full.clear();
    time_block.clear();
    time_overhead.clear();

  } // for
}
