// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeSmoothJacobi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputeSmoothJacobi class

#include "cello.hpp"

#include "enzo.hpp"

// #define DEBUG_SMOOTH
// #define DEBUG_VERBOSE

#if defined(DEBUG_SMOOTH) && defined(DEBUG_VERBOSE)
#   define PRINT_BDRX				\
  for (int ix=0; ix<mx; ix++) {			\
    for (int iy=0; iy<my; iy++) {		\
      for (int iz=0; iz<mz; iz++) {		\
	int i = ix + mx*(iy + my*iz);		\
	CkPrintf ("[%d %d] %d B=%6.3g D=%6.3g R=%6.3g X=%6.3g\n",	\
		ix,iy,i,B[i],D[i],R[i],X[i]);	\
      }						\
    }						\
  }
#else 
#  define PRINT_BDRX ;
#endif

//----------------------------------------------------------------------

EnzoComputeSmoothJacobi::EnzoComputeSmoothJacobi
( Matrix * A, int ix, int ib, int ir, int id, 
  double weight, int iter_max) throw()
  : A_ (A),
    ix_ (ix),
    ib_ (ib),
    ir_ (ir),
    id_ (id),
    w_(weight),
    n_(iter_max)
{
}

//----------------------------------------------------------------------

void EnzoComputeSmoothJacobi::compute ( Block * block) throw()
{
  Field field = block->data()->field();
  if (field.precision(0) == precision_single) {
    compute_<float>(block);
  } else if (field.precision(0) == precision_double) {
    compute_<double>(block);
  } else if (field.precision(0) == precision_quadruple) {
    compute_<long double>(block);
  }
}

//----------------------------------------------------------------------

template <typename T>
void EnzoComputeSmoothJacobi::compute_(Block * block)
{
  Field field = block->data()->field();

  int mx,my,mz;
  field.dimensions(ix_,&mx,&my,&mz);

  for (int iter=0; iter<n_; iter++) {
    
    /// Loop bounds minimal given iteration, ending at
    /// ghost_depth iter=n_-1

    //        const int g0 = MAX(gx - (n_-1-iter), 0);
    const int g0 = 1;

    const int ix0 = (mx > 1) ? g0 : 0;
    const int iy0 = (my > 1) ? g0 : 0;
    const int iz0 = (mz > 1) ? g0 : 0;

    // XXXX iter = 0
#ifdef DEBUG_SMOOTH
    printf ("%s:%d %s DEBUG_SMOOTH Computing residual %d\n",
	    __FILE__,__LINE__,block->name().c_str(),iter);
#endif
    PRINT_BDRX;

    A_->residual (ir_, ib_, ix_, block,g0);

#ifdef DEBUG_SMOOTH
    printf ("%s:%d %s DEBUG_SMOOTHComputing diagonal\n",
	    __FILE__,__LINE__,block->name().c_str());
#endif
    PRINT_BDRX;
    
    A_->diagonal (id_, block,g0);

#ifdef DEBUG_SMOOTH
    printf ("%s:%d %s DEBUG_SMOOTH Computing X=R/D\n",
	    __FILE__,__LINE__,block->name().c_str());
#endif
    PRINT_BDRX;

#ifdef DEBUG_SMOOTH
    printf ("%s:%d DEBUG LOOP LIMITS X=R/D ix=%d:%d iy=%d:%d iz=%d:%d\n",
	    __FILE__,__LINE__,ix0,mx-ix0,iy0,my-iy0, iz0, mz-iz0);
#endif

    T * D = (T*) field.values(id_);
    T * R = (T*) field.values(ir_);
    T * X = (T*) field.values(ix_);

    for (int iz=iz0; iz<mz-iz0; iz++) {
      for (int iy=iy0; iy<my-iy0; iy++) {
	for (int ix=ix0; ix<mx-ix0; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  X[i] += R[i] / D[i];
	}
      }
    }

#ifdef DEBUG_SMOOTH
    printf ("%s:%d %s DEBUG_SMOOTH  Done with SmoothJacobi\n",
	    __FILE__,__LINE__,block->name().c_str());
#endif
    PRINT_BDRX;

  }
}

