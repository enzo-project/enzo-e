// See LICENSE_CELLO file for license and copyright information

/// @file     field_ProlongLinear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of default linear prolongation

#include "problem.hpp"

//----------------------------------------------------------------------

ProlongLinear::ProlongLinear() throw()
  : Prolong ()
{
  TRACE("ProlongLinear::ProlongLinear");
}

//----------------------------------------------------------------------

void ProlongLinear::apply 
(
 precision_type precision,
 void * values_f,
 int ndx_f, int ndy_f, int ndf_z, 
 int nx_f,  int ny_f,  int nz_f,
 const void * values_c,
 int ndx_c, int ndy_c, int ndc_z, 
 int nx_c,  int ny_c,  int nz_c)
{
  switch (precision)  {
  case precision_single:
    apply_( (float *) values_f,
	    ndx_f, ndy_f, ndf_z, 
	    nx_f,  ny_f,  nz_f,
	    (const float *) values_c,
	    ndx_c, ndy_c, ndc_z, 
	    nx_c,  ny_c,  nz_c);
    break;
  case precision_double:
    apply_( (double *) values_f,
	    ndx_f, ndy_f, ndf_z, 
	    nx_f,  ny_f,  nz_f,
	    (const double *) values_c,
	    ndx_c, ndy_c, ndc_z, 
	    nx_c,  ny_c,  nz_c);
    break;

      }
}

template <class T>
void ProlongLinear::apply_
(
 T * values_f,
 int ndx_f, int ndy_f, int ndf_z, 
 int nx_f,  int ny_f,  int nz_f,
 const T * values_c,
 int ndx_c, int ndy_c, int ndc_z, 
 int nx_c,  int ny_c,  int nz_c)
{
  int dx_c = 1;
  int dy_c = ndx_c;
  int dz_c = ndy_c;

  const double c1[4] = { 5.0*0.25, 3.0*0.25, 1.0*0.25, -1.0*0.25};
  const double c2[4] = {-1.0*0.25, 1.0*0.25, 3.0*0.25,  5.0*0.25};

	  
  if (ny_f==1) {

    ASSERT ("ProlongLinear::apply_",
	    "fine array must be twice the size of the coarse array",
	    nx_f==nx_c*2);

    ASSERT ("ProlongLinear::apply_",
	    "fine grid array sizes must be divisible by 4",
	    nx_f % 4 == 0);

    for (int ix0=0; ix0<nx_f; ix0+=4) {
      int i_c = ix0/2;
      for (int ix=ix0; ix<ix0+4; ix++) {
	int icx = ix-ix0;
	int i_f = ix;
	values_f[i_f] = ( c1[icx]*values_c[i_c] + 
			  c2[icx]*values_c[i_c+dx_c]);
      }
    }

  } else if (nz_f == 1) {

    ASSERT ("ProlongLinear::apply_",
	    "fine array must be twice the size of the coarse array",
	    ny_f==ny_c*2);

    ASSERT ("ProlongLinear::apply_",
	    "fine grid array sizes must be divisible by 4",
	    ny_f % 4 == 0);

    for (int ix0=0; ix0<nx_f; ix0+=4) {
      for (int iy0=0; iy0<ny_f; iy0+=4) {

	int i_c = ix0/2 + ndy_c*(iy0/2);

	for (int ix=ix0; ix<ix0+4; ix++) {
	  int icx = ix-ix0;
	  for (int iy=iy0; iy<iy0+4; iy++) {
	    int icy = iy-iy0;

	    int i_f = ix + ndx_f*iy;

	    values_f[i_f] = 
	      ( c1[icx]*c1[icy]*values_c[i_c] +
		c2[icx]*c1[icy]*values_c[i_c+dx_c] +
		c1[icx]*c2[icy]*values_c[i_c     +dy_c] +
		c2[icx]*c2[icy]*values_c[i_c+dx_c+dy_c]);
	  }
	}
      }
    }
  } else {

    ASSERT ("ProlongLinear::apply_",
	    "fine array must be twice the size of the coarse array",
	    nz_f==nz_c*2);

    ASSERT ("ProlongLinear::apply_",
	    "fine grid array sizes must be divisible by 4",
	    nz_f % 4 == 0);

    for (int ix0=0; ix0<nx_f; ix0+=4) {
      for (int iy0=0; iy0<ny_f; iy0+=4) {
	for (int iz0=0; iz0<nz_f; iz0+=4) {

	  int i_c = ix0/2 + ndx_c*(iy0/2 + ndy_c*iz0/2);

	  for (int ix=ix0; ix<ix0+4; ix++) {
	    int icx = ix-ix0;
	    for (int iy=iy0; iy<iy0+4; iy++) {
	      int icy = iy-iy0;
	      for (int iz=iz0; iz<iz0+4; iz++) {
		int icz = iz-iz0;

		int i_f = ix + ndx_f*(iy + ndy_f*iz);

		values_f[i_f] = 
		  ( c1[icx]*c1[icy]*c1[icz]*values_c[i_c] +
		    c2[icx]*c1[icy]*c1[icz]*values_c[i_c+dx_c] +
		    c1[icx]*c2[icy]*c1[icz]*values_c[i_c     +dy_c] +
		    c2[icx]*c2[icy]*c1[icz]*values_c[i_c+dx_c+dy_c] +
		    c1[icx]*c1[icy]*c2[icz]*values_c[i_c          +dz_c] +
		    c2[icx]*c1[icy]*c2[icz]*values_c[i_c+dx_c     +dz_c] +
		    c1[icx]*c2[icy]*c2[icz]*values_c[i_c     +dy_c+dz_c] +
		    c2[icx]*c2[icy]*c2[icz]*values_c[i_c+dx_c+dy_c+dz_c]);
	      }
	    }
	  }
	}
      }
    }
  }
}

//======================================================================

