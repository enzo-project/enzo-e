//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      time_matvec.cpp
 * @brief     Program timing sparse matrix-vector multiplies
 * @author    James Bordner
 * @date      Thu Apr 24 16:25:43 PDT 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "cello.h"
#include "error.hpp"
#include "performance.hpp"

int main(int argc, char ** argv)
{
  Timer time_overhead;
  Timer time_const;
  Timer time_symm;
  Timer time_full;
  Timer time_block;
  int i,j,k;
  int i0,i1,i2;

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

  printf ("size  mflop-const mflop-symm mflop-full mflop-block\n");

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
    Scalar ac[7] = {1,1,1,1,1,1,1};

    {
      if (a0 && axp && axm && ayp && aym && azp && azm) {

	for (j=0; j<count; j++) {

	  time_const.start();

	  for (i0=1; i0<k+1; i0++) {
	    for (i1=1; i1<k+1; i1++) {
	      for (i2=1; i2<k+1; i2++) {
		i = i2 + k2*(i1 + k2*i0);
		b[i] = ac[0]*x[i-k22]
		  +    ac[1]*x[i-k2]
		  +    ac[2]*x[i-1]
		  +    ac[3]*x[i]
		  +    ac[4]*x[i+1]
		  +    ac[5]*x[i+k2]
		  +    ac[6]*x[i+k22];
	      }
	    }
	  }

	  time_const.stop();

	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 0;
	  }
	}

	for (j=0; j<count; j++) {
	  time_overhead.start();
	  time_overhead.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 0;
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
		i = i2 + k2*(i1 + k2*i0);
		b[i] = azp[i-k22]*x[i-k22]
		  +    ayp[i-k2] *x[i-k2]
		  +    axp[i-1]  *x[i-1]
		  +     a0[i]      *x[i]
		  +    axp[i]    *x[i+1]
		  +    ayp[i]   *x[i+k2]
		  +    azp[i]*x[i+k22];
	      }
	    }
	  }
	  time_symm.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 0;
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
		i = i2 + k2*(i1 + k2*i0);
		b[i] = azm[i]*x[i-k22]
		  +    aym[i]*x[i-k2]
		  +    axm[i]*x[i-1]
		  +     a0[i]*x[i]
		  +    axp[i]*x[i+1]
		  +    ayp[i]*x[i+k2]
		  +    azp[i]*x[i+k22];
	      }
	    }
	  }
	  time_full.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 0;
	  }
	}
      }
    }

    {
      if (a0 && axp && axm && ayp && aym && azp && azm) {

	// Shift b to improve inter-array cache effects
	// Force k2 to be odd to improve intra-array cache effects.

	int adjust = (1 - k2%2);

	b   = b  + adjust*((k2+1)*(k2+1)*(k2+1) - (k2*k2*k2));
	k2  = k2 + adjust;
	k22 = k2*k2;


	for (j=0; j<count; j++) {

	  time_block.start();

	  for (i0=1; i0<k+1; i0++) {
	    for (i1=1; i1<k+1; i1++) {
	      for (i2=1; i2<k+1; i2++) {
		i = i2 + k2*(i1 + k2*i0);
		b[i] = ac[0]*x[i-k22]
		  +    ac[1]*x[i-k2]
		  +    ac[2]*x[i-1]
		  +    ac[3]*x[i]
		  +    ac[4]*x[i+1]
		  +    ac[5]*x[i+k2]
		  +    ac[6]*x[i+k22];
	      }
	    }
	  }

	  time_block.stop();

	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 0;
	  }
	}

	for (j=0; j<count; j++) {
	  time_overhead.start();
	  time_overhead.stop();
	  for (i=0; i<CACHE_SIZE; i++) {
	    cache[i] = 0;
	  }
	}
      }
    }

    int num_flops = count*k*k*k*13+int(cache[0]);
    printf ("%d %10.2f %10.2f %10.2f %10.2f\n",
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
