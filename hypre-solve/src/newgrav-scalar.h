//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Scalar type definitions

 /**
 *********************************************************************
 *
 * @file      scalar.hpp
 * @brief     Scalar related defines
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 *********************************************************************
 */

typedef double Scalar ;

#define SCALAR_MAX 1e38;
#define SCALAR_SCANF "%lf"
#define SCALAR_PRINTF "%le "
#define MPI_SCALAR MPI_DOUBLE

//----------------------------------------------------------------------
// TEMPORARY
//----------------------------------------------------------------------

#define  WRITE_B_SUM(grid) \
  { \
    int n0,n1,n2; \
    Scalar * f; \
    f = grid->get_f(&n0,&n1,&n2);			\
    double sum = 0.0; \
    for (int i2=0; i2<n2; i2++) { \
      for (int i1=0; i1<n1; i1++) { \
	for (int i0=0; i0<n0; i0++) { \
	  int k = i0 + n0*(i1 + n1*i2); \
	  sum += f[k]; \
	} \
      } \
    } \
    printf ("%s:%d %d DEBUG sum(%d %d %d)=%g\n",__FILE__,__LINE__,pmpi->ip(),n0,n1,n2,sum); \
  }
//----------------------------------------------------------------------



