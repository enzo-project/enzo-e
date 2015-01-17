// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlongPoisson.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-01-16
/// 
/// @brief    [\ref Enzo] Implentation of "Poisson" prolongation.
///
/// The EnzoProlongPoisson class is used for interpolation of coarse
/// Field values to neighboring fine Block's in the
/// EnzoMethodGravityCg method.  It is used to ensure that the first
/// derivative as well as values of the interpolated functions are
/// continuous.  It does this by the following:
///
///    - 1D: coarse value is sent directly, which on the receiving
///      side is used to perform a quadratic interpolation to find the
///      ghost value
///    - 2D: three coarse values are quadratically interpolated to be
///      in line with fine ghost zones and send to neighbor's ghost
///      zones.  At receiving end quadratic interpolation is used to
///      find the ghost value as in the 1D case.
///    - 3D: three coarse values are bilinearly interpolated to be in
///      line with fine ghost zones and send to neighbor's ghost
///      zones.  At receiving end quadratic interpolation is used to
///      find the ghost value as in the 1D case.

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoProlongPoisson::EnzoProlongPoisson() throw()
  : Prolong ()
{
  printf ("EnzoProlongPoisson\n");
  TRACE("EnzoProlongPoisson::EnzoProlongPoisson");
}

//----------------------------------------------------------------------

int EnzoProlongPoisson::apply 
( precision_type precision,
  void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
  const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3])
{
  TRACE6("EnzoProlongPoisson fine   %d:%d %d:%d %d:%d",
	 im3_f[0],n3_f[0]+im3_f[0],
	 im3_f[1],n3_f[1]+im3_f[1],
	 im3_f[2],n3_f[2]+im3_f[2]);

  TRACE6("EnzoProlongPoisson coarse %d:%d %d:%d %d:%d",
	 im3_c[0],n3_c[0]+im3_c[0],
	 im3_c[1],n3_c[1]+im3_c[1],
	 im3_c[2],n3_c[2]+im3_c[2]);

  switch (precision)  {

  case precision_single:

    return apply_(       (float *) values_f, nd3_f, im3_f, n3_f,
		   (const float *) values_c, nd3_c, im3_c, n3_c);

    break;

  case precision_double:

    return apply_(       (double *) values_f, nd3_f, im3_f, n3_f,
		   (const double *) values_c, nd3_c, im3_c, n3_c);

    break;

  default:

    ERROR1 ("EnzoProlongPoisson::apply()",
	    "Unknown precision %d",
	    precision);

    return 0;
  }
}

//----------------------------------------------------------------------

template <class T>
int EnzoProlongPoisson::apply_
(       T * values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
  const T * values_c, int nd3_c[3], int im3_c[3], int n3_c[3])
{
  int dx_c = 1;
  int dy_c = nd3_c[0];
  int dz_c = nd3_c[1];

  const double c1[4] = { 5.0*0.25, 3.0*0.25, 1.0*0.25, -1.0*0.25};
  const double c2[4] = {-1.0*0.25, 1.0*0.25, 3.0*0.25,  5.0*0.25};

  int rank = (nd3_f[2] > 1) ? 3 : ( (nd3_f[1] > 1) ? 2 : 1 );

  for (int i=0; i<rank; i++) {
    const char * xyz = "xyz";
    ASSERT3 ("EnzoProlongPoisson::apply_",
	     "fine array %c-axis %d must be twice the size of the coarse axis %d",
	     xyz[i],n3_c[i],n3_f[i],
	     n3_f[i]==n3_c[i]*2);

    ASSERT2 ("EnzoProlongPoisson::apply_",
	     "fine grid %c-axis %d must be divisible by 4",
	     xyz[i],n3_f[i],n3_f[i] % 4 == 0);
  }

  if (n3_f[1]==1) {

    for (int ix0=0; ix0<n3_f[0]; ix0+=4) {
      int ic_x = ix0/2;
      for (int ix=ix0; ix<ix0+4; ix++) {
	int icx = ix-ix0;
	int if_x = ix;

	int i_c = im3_c[0]+ic_x;
	int i_f = im3_f[0]+if_x;

	values_f[i_f] = ( c1[icx]*values_c[i_c] + 
			  c2[icx]*values_c[i_c+dx_c]);

	if (positive_ && (values_f[i_f] < 0)) {
	  // revert to poisson
	  int icc = i_c + (icx/2)*dx_c;
	  values_f[i_f] = values_c[icc];	  
	}
      }
    }

    return (sizeof(T) * n3_c[0]);


  } else if (n3_f[2] == 1) {

    for (int ix0=0; ix0<n3_f[0]; ix0+=4) {
      int ic_x = ix0/2;
      for (int iy0=0; iy0<n3_f[1]; iy0+=4) {
	int ic_y = iy0/2;
	for (int ix=ix0; ix<ix0+4; ix++) {
	  int icx = ix-ix0;
	  int if_x = ix;
	  for (int iy=iy0; iy<iy0+4; iy++) {
	    int icy = iy-iy0;
	    int if_y = iy;

	    int i_c = (im3_c[0]+ic_x) + nd3_c[0]*
	      (       (im3_c[1]+ic_y));
	    int i_f = (im3_f[0]+if_x) + nd3_f[0]*
	      (       (im3_f[1]+if_y));

	    values_f[i_f] = 
	      ( c1[icx]*c1[icy]*values_c[i_c] +
		c2[icx]*c1[icy]*values_c[i_c+dx_c] +
		c1[icx]*c2[icy]*values_c[i_c     +dy_c] +
		c2[icx]*c2[icy]*values_c[i_c+dx_c+dy_c]);
	    if (positive_ && (values_f[i_f] < 0)) {
	      // revert to poisson
	      int icc = i_c + (icx/2)*dx_c + (icy/2)*dy_c;
	      values_f[i_f] = values_c[icc];	  
	    }
	  }
	}
      }
    }
    return (sizeof(T) * n3_c[0]*n3_c[1]);

  } else {

    for (int ix0=0; ix0<n3_f[0]; ix0+=4) {
      int ic_x = ix0/2;
      for (int iy0=0; iy0<n3_f[1]; iy0+=4) {
	int ic_y = iy0/2;
	for (int iz0=0; iz0<n3_f[2]; iz0+=4) {
	  int ic_z = iz0/2;

	  for (int ix=ix0; ix<ix0+4; ix++) {
	    int icx = ix-ix0;
	    int if_x = ix;
	    for (int iy=iy0; iy<iy0+4; iy++) {
	      int icy = iy-iy0;
	      int if_y = iy;
	      for (int iz=iz0; iz<iz0+4; iz++) {
		int icz = iz-iz0;
		int if_z = iz;

		int i_c = (im3_c[0]+ic_x) + nd3_c[0]*
		  (       (im3_c[1]+ic_y) + nd3_c[1]*
			  (im3_c[2]+ic_z));
		int i_f = (im3_f[0]+if_x) + nd3_f[0]*
		  (       (im3_f[1]+if_y) + nd3_f[1]*
			  (im3_f[2]+if_z));

		values_f[i_f] = 
		  ( c1[icx]*c1[icy]*c1[icz]*values_c[i_c] +
		    c2[icx]*c1[icy]*c1[icz]*values_c[i_c+dx_c] +
		    c1[icx]*c2[icy]*c1[icz]*values_c[i_c     +dy_c] +
		    c2[icx]*c2[icy]*c1[icz]*values_c[i_c+dx_c+dy_c] +
		    c1[icx]*c1[icy]*c2[icz]*values_c[i_c          +dz_c] +
		    c2[icx]*c1[icy]*c2[icz]*values_c[i_c+dx_c     +dz_c] +
		    c1[icx]*c2[icy]*c2[icz]*values_c[i_c     +dy_c+dz_c] +
		    c2[icx]*c2[icy]*c2[icz]*values_c[i_c+dx_c+dy_c+dz_c]);
		if (positive_ && (values_f[i_f] < 0)) {
		  // revert to poisson
		  int icc = i_c + (icx/2)*dx_c + (icy/2)*dy_c + (icz/2)*dz_c;
		  values_f[i_f] = values_c[icc];	  
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return (sizeof(T) * n3_c[0]*n3_c[1]*n3_c[2]);

}

//======================================================================

