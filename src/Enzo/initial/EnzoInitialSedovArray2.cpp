// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialSedovArray2.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Sat Jun 15 13:44:27 PDT 2013
/// @brief    Implementation of an array of Sedov problems, one per Block
///
/// This problem is designed for weak-scaling studies.  Each block
/// contains a small Sedov blast problem.  Problem-specific parameters
/// are used to define the initial values, and whether all blasts go
/// at once or intermitently

#include "Enzo/initial/initial.hpp"
#include "Enzo/enzo.hpp"
#include "Cello/cello.hpp"

//----------------------------------------------------------------------

EnzoInitialSedovArray2::EnzoInitialSedovArray2 
(const EnzoConfig * config) throw ()
: Initial(config->initial_cycle, config->initial_time)
{
  array_[0]        = config->initial_sedov_array[0];
  array_[1]        = config->initial_sedov_array[1];
  radius_relative_ = config->initial_sedov_radius_relative;
  pressure_in_     = config->initial_sedov_pressure_in;
  pressure_out_    = config->initial_sedov_pressure_out;
  density_         = config->initial_sedov_density;

}

//----------------------------------------------------------------------

void EnzoInitialSedovArray2::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  Initial::pup(p);

  TRACEPUP;

  PUParray(p,array_,2);
  p | radius_relative_;
  p | pressure_in_;
  p | pressure_out_;
  p | density_;
}

//----------------------------------------------------------------------
void EnzoInitialSedovArray2::enforce_block
( Block * block, const Hierarchy  * hierarchy ) throw()

{

  if (!block->is_leaf()) {
    block->initial_done();
    return;
  }

  ASSERT("EnzoInitialSedovArray2",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  ASSERT("EnzoInitialSedovArray2",
	 "Insufficient number of fields",
	 field.field_count() >= 4);

  enzo_float *  d = (enzo_float *) field.values ("density");
  enzo_float * te = (enzo_float *) field.values ("total_energy");

  int nx,ny;
  field.size(&nx,&ny);

  double xbm,ybm;
  block->data()->lower(&xbm,&ybm);

  double xbp,ybp;
  block->data()->upper(&xbp,&ybp);

  double hx,hy;
  field.cell_width(xbm,xbp,&hx,
			  ybm,ybp,&hy);

  // Parameters

  const double sedov_radius = radius_relative_/array_[0];
  const double sedov_radius_2 = sedov_radius*sedov_radius;

  const enzo_float gamma = enzo::fluid_props()->gamma();

  const double sedov_te_in = pressure_in_  / ((gamma - 1.0) * density_);
  const double sedov_te_out= pressure_out_ / ((gamma - 1.0) * density_);

  int gx,gy;
  field.ghost_depth(0,&gx,&gy);

  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;

  // clear all fields

  for (int iv=0; iv<field.field_count(); iv++) {

    enzo_float * array = (enzo_float *) field.values (iv);

    for (int iy=0; iy<ndy; iy++) {
      for (int ix=0; ix<ndx; ix++) {
	int i = INDEX2(ix,iy,ndx);
	array[i] = 0.0;
      }
    }
  }

  // background 

  for (int iy=0; iy<ndy; iy++) {
    for (int ix=0; ix<ndx; ix++) {

      int i = INDEX2(ix,iy,ndx);
      d[i]  = density_;
      te[i] = sedov_te_out;

    }
  }

  // array of explosions

  // (kx,ky) index of explosion in domain
  // (x,y) position in block
  int kxm = 0;
  int kym = 0;
  int kxp = array_[0];
  int kyp = array_[1];
  TRACE2 ("SEDOV: %d %d",kxp,kyp);

  double xdm,ydm;
  hierarchy->lower(&xdm,&ydm);
  double xdp,ydp;
  hierarchy->upper(&xdp,&ydp);

  double hxa = (xdp-xdm) / array_[0];
  double hya = (ydp-ydm) / array_[1];

  for (int ky=kym; ky<kyp; ky++) {
    double yc = hya*(0.5+ky);
    for (int kx=kxm; kx<kxp; kx++) {
      double xc = hxa*(0.5+kx);
      TRACE2("xc,yc = %f %f %f",xc,yc);

      // (explosion center xc,yc)

      for (int iy=0; iy<ndy; iy++) {
	double y = ybm + (iy - gy + 0.5)*hy - yc;
	for (int ix=0; ix<ndx; ix++) {
	  double x = xbm + (ix - gx + 0.5)*hx - xc;
	  double r2 = x*x + y*y;

	  int i = INDEX2(ix,iy,ndx);
	  
	  if (r2 < sedov_radius_2) {
	    te[i] = sedov_te_in;
	  }
	}
      }
    }
  }

  block->initial_done();
}

