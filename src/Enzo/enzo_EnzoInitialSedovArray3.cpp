// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialSedovArray3.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of an array of Sedov problems, one per Block
///
/// This problem is designed for weak-scaling studies.  Each block
/// contains a small Sedov blast problem.  Problem-specific parameters
/// are used to define the initial values, and whether all blasts go
/// at once or intermitently

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialSedovArray3::EnzoInitialSedovArray3 
(const EnzoConfig * config) throw ()
: Initial(config->initial_cycle, config->initial_time) 
{
  array_[0] = config->sedov_array[0];
  array_[1] = config->sedov_array[1];
  array_[2] = config->sedov_array[2];
  radius_relative_ = config->sedov_radius_relative;
  pressure_in_     = config->sedov_pressure_in;
  pressure_out_    = config->sedov_pressure_out;
  density_         = config->sedov_density;
}

//----------------------------------------------------------------------

void EnzoInitialSedovArray3::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,array_,3);
  p | radius_relative_;
  p | pressure_in_;
  p | pressure_out_;
  p | density_;
  
}

//----------------------------------------------------------------------
void EnzoInitialSedovArray3::enforce_block
(
 CommBlock * comm_block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  ASSERT("EnzoInitialSedovArray3",
	 "CommBlock does not exist",
	 comm_block != NULL);

  FieldBlock * field_block = comm_block->block()->field_block();

  ASSERT("EnzoInitialSedovArray3",
	 "Insufficient number of fields",
	 field_descr->field_count() >= 4);

  enzo_float *  d = (enzo_float *) field_block->field_values
    (field_descr->field_id("density"));
  
  enzo_float * te = 0;

  te = (enzo_float *) field_block->field_values
    (field_descr->field_id("total_energy"));

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  double xbm,ybm,zbm;
  comm_block->block()->lower(&xbm,&ybm,&zbm);

  double xbp,ybp,zbp;
  comm_block->block()->upper(&xbp,&ybp,&zbp);

  double hx,hy,hz;
  field_block->cell_width(xbm,xbp,&hx,
			  ybm,ybp,&hy,
			  zbm,zbp,&hz);

  // Parameters

  const double sedov_radius = radius_relative_/array_[0];
  const double sedov_radius_2 = sedov_radius*sedov_radius;

  const double sedov_te_in = 
    pressure_in_  / ((EnzoBlock::Gamma - 1.0) * density_);
  const double sedov_te_out= 
    pressure_out_ / ((EnzoBlock::Gamma - 1.0) * density_);

  int gx,gy,gz;
  field_descr->ghosts(0,&gx,&gy,&gz);

  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;
  int ndz = nz + 2*gz;

  // clear all fields

  for (int iv=0; iv<field_descr->field_count(); iv++) {

    enzo_float * field = (enzo_float *) field_block->field_values (iv);

    for (int iz=0; iz<ndz; iz++) {
      for (int iy=0; iy<ndy; iy++) {
	for (int ix=0; ix<ndx; ix++) {
	  int i = INDEX(ix,iy,iz,ndx,ndy);
	  field[i] = 0.0;
	}
      }
    }
  }

  // background 

  for (int iz=0; iz<ndz; iz++) {
    for (int iy=0; iy<ndy; iy++) {
      for (int ix=0; ix<ndx; ix++) {

	int i = INDEX(ix,iy,iz,ndx,ndy);
	d[i]  = density_;
	te[i] = sedov_te_out;

      }
    }
  }

  // array of explosions

  // (kx,ky,kz) index of explosion in domain
  // (x,y,z) position in block
  int kxm = 0;
  int kym = 0;
  int kzm = 0;
  int kxp = array_[0];
  int kyp = array_[1];
  int kzp = array_[2];
  TRACE3 ("SEDOV: %d %d %d",kxp,kyp,kzp);

  double xdm,ydm,zdm;
  hierarchy->lower(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp;
  hierarchy->upper(&xdp,&ydp,&zdp);

  double hxa = (xdp-xdm) / array_[0];
  double hya = (ydp-ydm) / array_[1];
  double hza = (zdp-zdm) / array_[2];

  for (int kz=kzm; kz<kzp; kz++) {
    double zc = hza*(0.5+kz);
    for (int ky=kym; ky<kyp; ky++) {
      double yc = hya*(0.5+ky);
      for (int kx=kxm; kx<kxp; kx++) {
	double xc = hxa*(0.5+kx);
	TRACE3("xc,yc,zc = %f %f %f",xc,yc,zc);

	// (explosion center xc,yc,zc)

	for (int iz=0; iz<ndz; iz++) {
	  double z = zbm + (iz - gz + 0.5)*hz - zc;
	  for (int iy=0; iy<ndy; iy++) {
	    double y = ybm + (iy - gy + 0.5)*hy - yc;
	    for (int ix=0; ix<ndx; ix++) {
	      double x = xbm + (ix - gx + 0.5)*hx - xc;
	      double r2 = x*x + y*y + z*z;

	      int i = INDEX(ix,iy,iz,ndx,ndy);
	  
	      if (r2 < sedov_radius_2) {
		te[i] = sedov_te_in;
	      }
	    }
	  }
	}
      }
    }
  }
  
}

