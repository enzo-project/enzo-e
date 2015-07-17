// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Matrix.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include "compute.hpp"

// #define DEBUG_MATRIX

//----------------------------------------------------------------------

void Matrix::residual (int ir, int ib, int ix, Block * block, int g0) throw()
{

  matvec(ir,ix,block);

  Field field = block->data()->field();

  void * B = field.values(ib);
  void * R = field.values(ir);

  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);

  int precision = field.precision(0);

  if      (precision == precision_single)    
    residual_((float *)(R), (float *)(B),
	      mx,my,mz,g0);
  else if (precision == precision_double)    
    residual_((double *)(R), (double *)(B),
	      mx,my,mz,g0);
  else if (precision == precision_quadruple) 
    residual_((long double *)(R), (long double *)(B),
	      mx,my,mz,g0);
  else 
    ERROR1("Matrix::residual()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void Matrix::residual_ (T * r, T * b,
			int mx, int my, int mz,
			int g0) throw()
{

  const int ix0 = (mx > 1) ? g0 : 0;
  const int iy0 = (my > 1) ? g0 : 0;
  const int iz0 = (mz > 1) ? g0 : 0;

#ifdef DEBUG_MATRIX
    printf ("%s:%d DEBUG_LOOP_LIMITS: (%d:%d) (%d:%d) (%d:%d)\n",
	    __FILE__,__LINE__, ix0,mx-ix0,iy0,my-iy0,iz0,mz-iz0);
#endif

  for (int iz=iz0; iz<mz-iz0; iz++) {
    for (int iy=iy0; iy<my-iy0; iy++) {
      for (int ix=ix0; ix<mx-ix0; ix++) {

	const int i=ix + mx*(iy + mz*iz);

	r[i] = b[i] - r[i];
      }
    }
  }
}
//======================================================================

