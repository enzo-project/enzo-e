// See LICENSE_CELLO file for license and copyright information

/// @file     field_RestrictLinear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    Implentation of default linear restriction

#include "problem.hpp"

//----------------------------------------------------------------------

RestrictLinear::RestrictLinear() throw()
  : Restrict ()
{
  TRACE("RestrictLinear::RestrictLinear");
}

//----------------------------------------------------------------------

void RestrictLinear::apply 
(
 FieldBlock        * field_block_c, 
 const  FieldBlock * field_block_f, 
 const FieldDescr * field_descr,
 int icx, int icy, int icz)
{
  INCOMPLETE("RestructLinear::apply");

  return;

  for (int index=0; index<field_descr->field_count(); index++) {

    int gx,gy,gz;
    field_descr->ghosts(index,&gx,&gy,&gz);

    int nx,ny,nz;
    field_block_f->size(&nx,&ny,&nz);

    int ndx = nx + 2*gx;
    int ndy = ny + 2*gy;
    int ndz = nz + 2*gz;

    int ixm = icx * nx/2;
    int iym = icy * ny/2;
    int izm = icz * nz/2;

    void *       values_c = field_block_c->field_values(index);
    const void * values_f = field_block_f->field_values(index);

    
    switch (field_descr->precision(index)) {

    case precision_single:
      
      interpolate_((float *) values_c,
		   (const float *) values_f,
		   ndx,ndy,ndz,
		   ixm,iym,izm,
		   nx,ny,nz,
		   gx,gy,gz);
      break;

    case precision_double:
      
      interpolate_((double *) values_c,
		   (const double *) values_f,
		   ndx,ndy,ndz,
		   ixm,iym,izm,
		   nx, ny, nz,
		   gx,gy,gz);
      break;

    default:

      ERROR2 ("RestrictLinear::apply()",
	      "Unknown precision %d for field %d",
	      field_descr->precision(index),index);

    }
  }
}

//----------------------------------------------------------------------

template<class T>
void RestrictLinear::interpolate_(T * values_c,
				 const T * values_f,
				 int ndx, int ndy, int ndz,
				 int ixm, int iym, int izm,
				 int nx, int ny, int nz,
				 int gx, int gy, int gz)
{
  INCOMPLETE("RestrictLinear::interpolate_()");
  return;
  int dx = 1;
  int dy = ndx;
  int dz = ndy;
  double c1[4] = { 5.0*0.25, 3.0*0.25, 1.0*0.25, -1.0*0.25};
  double c2[4] = {-1.0*0.25, 1.0*0.25, 3.0*0.25,  5.0*0.25};
  if (ny==1) {
    ASSERT ("RestrictLinear::interpolate_",
	    "Block sizes must be divisible by 4",
	    nx % 4 == 0);

    for (int ix0=0; ix0<nx; ix0+=4) {
      int iH = ixm + ix0/2;
      for (int ix=ix0; ix<ix0+4; ix++) {
	int icx = ix-ix0;
	int ih = ix;
	values_c[ih] = ( c1[icx]*values_c[iH] + c2[icx]*values_c[iH+dx]);
      }
    }

  } else if (nz == 1) {
    for (int ix0=0; ix0<nx; ix0+=4) {
      for (int iy0=0; iy0<ny; iy0+=4) {

	int iH = ixm + ix0/2 + ndy*(iym + iy0/2);

	for (int ix=ix0; ix<ix0+4; ix++) {
	  int icx = ix-ix0;
	  for (int iy=iy0; iy<iy0+4; iy++) {
	    int icy = iy-iy0;

	    int ih = ix + ndx*iy;

	    values_c[ih] = 
	      ( c1[icx]*c1[icy]*values_c[iH] +
		c2[icx]*c1[icy]*values_c[iH+dx] +
		c1[icx]*c2[icy]*values_c[iH+dy] +
		c2[icx]*c2[icy]*values_c[iH+dy+dy]);
	  }
	}
      }
    }
  } else {

    for (int ix0=0; ix0<nx; ix0+=4) {
      for (int iy0=0; iy0<ny; iy0+=4) {
	for (int iz0=0; iz0<nz; iz0+=4) {

	  int iH = ixm + ix0/2 + ndy*(iym + iy0/2);

	  for (int ix=ix0; ix<ix0+4; ix++) {
	    int icx = ix-ix0;
	    for (int iy=iy0; iy<iy0+4; iy++) {
	      int icy = iy-iy0;
	      for (int iz=iz0; iz<iz0+4; iz++) {
		int icz = iz-iz0;

		int ih = ix + ndx*(iy + ndy*iz);

		values_c[ih] = 
		  ( c1[icx]*c1[icy]*c1[icz]*values_c[iH] +
		    c2[icx]*c1[icy]*c1[icz]*values_c[iH+dx] +
		    c1[icx]*c2[icy]*c1[icz]*values_c[iH+dy] +
		    c2[icx]*c2[icy]*c1[icz]*values_c[iH+dy+dy] +
		    c1[icx]*c1[icy]*c2[icz]*values_c[iH+dz] +
		    c2[icx]*c1[icy]*c2[icz]*values_c[iH+dx+dz] +
		    c1[icx]*c2[icy]*c2[icz]*values_c[iH+dy+dz] +
		    c2[icx]*c2[icy]*c2[icz]*values_c[iH+dy+dy+dz]);
	      }
	    }
	  }
	}
      }
    }
  }
}

//======================================================================

