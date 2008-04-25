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

main()
{
  Timer timer;
  int k,j,i0,i1,i2;

  for (k=8; k<=192; k++) {

    int k2 = k + 2;
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
    timer.start();

    int count = 10000000/k/k/k+1;
    for (j=0; j<count; j++) {
    for (i0=1; i0<k+1; i0++) {
      for (i1=1; i1<k+1; i1++) {
	for (i2=1; i2<k+1; i2++) {
	  int i = i0 + k2*(i1 + k2*i2);
	  b[i] = azm[i-k2*k2]*x[i-k2*k2]
	    +    aym[i-k2]   *x[i-k2]
	    +    axm[i-1]    *x[i-1]
	    +     a0[i]      *x[i]
	    +    axm[i+1]    *x[i+1]
	    +    aym[i+k2]   *x[i+k2]
	    +    azm[i+k2*k2]*x[i+k2*k2];
	}
      }
    }
    }

    timer.stop();
    int num_flops = count*k*k*k*13;
    printf ("%d %20.12f %20.12f \n",k,timer.value()/count, 10e-6*num_flops/timer.value()); 
    fflush(stdout);
    timer.clear();
    } else {
    printf ("%d 0\n",k); 
    fflush(stdout);
    }

  }
}
