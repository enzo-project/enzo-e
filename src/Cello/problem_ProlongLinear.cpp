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
 CommBlock        * comm_block_h, 
 const  CommBlock * comm_block_H, 
 const FieldDescr * field_descr,
 int icx, int icy, int icz)
{
  Block * block_h = comm_block_h->block();
  const Block * block_H = comm_block_H->block();
  FieldBlock * field_block_h = block_h->field_block();
  const FieldBlock * field_block_H = block_H->field_block();

  for (int index=0; index<field_descr->field_count(); index++) {

    int gx,gy,gz;
    field_descr->ghosts(index,&gx,&gy,&gz);

    int nx,ny,nz;
    field_block_h->size(&nx,&ny,&nz);

    int ndx = nx + 2*gx;
    int ndy = ny + 2*gy;
    int ndz = nz + 2*gz;

    int ixm = icx * nx/2;
    int iym = icy * ny/2;
    int izm = icz * nz/2;

    void * values_h       = field_block_h->field_values(index);
    const void * values_H = field_block_H->field_values(index);

    
    switch (field_descr->precision(index)) {

    case precision_single:
      
      interpolate_((float *) values_h,
		   (const float *) values_H,
		   ndx,ndy,ndz,
		   ixm,iym,izm,
		   nx,ny,nz,
		   gx,gy,gz);
      break;

    case precision_double:
      
      interpolate_((double *) values_h,
		   (const double *) values_H,
		   ndx,ndy,ndz,
		   ixm,iym,izm,
		   nx, ny, nz,
		   gx,gy,gz);
      break;

    default:

      ERROR2 ("ProlongLinear::apply()",
	      "Unknown precision %d for field %d",
	      field_descr->precision(index),index);

    }
  }
}

//----------------------------------------------------------------------

template<class T>
void ProlongLinear::interpolate_(T * values_h,
				 const T * values_H,
				 int ndx, int ndy, int ndz,
				 int ixm, int iym, int izm,
				 int nx, int ny, int nz,
				 int gx, int gy, int gz)
{

  int dx = 1;
  int dy = ndx;
  int dz = ndy;
  double c1[4] = { 5.0*0.25, 3.0*0.25, 1.0*0.25, -1.0*0.25};
  double c2[4] = {-1.0*0.25, 1.0*0.25, 3.0*0.25,  5.0*0.25};
  if (ny==1) {
    ASSERT ("ProlongLinear::interpolate_",
	    "Block sizes must be divisible by 4",
	    nx % 4 == 0);

    for (int ix0=0; ix0<nx; ix0+=4) {
      int iH = ixm + ix0/2;
      for (int ix=ix0; ix<ix0+4; ix++) {
	int icx = ix-ix0;
	int ih = ix;
	values_h[ih] = ( c1[icx]*values_H[iH] + c2[icx]*values_H[iH+dx]);
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

	    values_h[ih] = 
	      ( c1[icx]*c1[icy]*values_H[iH] +
		c2[icx]*c1[icy]*values_H[iH+dx] +
		c1[icx]*c2[icy]*values_H[iH+dy] +
		c2[icx]*c2[icy]*values_H[iH+dy+dy]);
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

		values_h[ih] = 
		  ( c1[icx]*c1[icy]*c1[icz]*values_H[iH] +
		    c2[icx]*c1[icy]*c1[icz]*values_H[iH+dx] +
		    c1[icx]*c2[icy]*c1[icz]*values_H[iH+dy] +
		    c2[icx]*c2[icy]*c1[icz]*values_H[iH+dy+dy] +
		    c1[icx]*c1[icy]*c2[icz]*values_H[iH+dz] +
		    c2[icx]*c1[icy]*c2[icz]*values_H[iH+dx+dz] +
		    c1[icx]*c2[icy]*c2[icz]*values_H[iH+dy+dz] +
		    c2[icx]*c2[icy]*c2[icz]*values_H[iH+dy+dy+dz]);
	      }
	    }
	  }
	}
      }
    }
  }
}

//======================================================================

