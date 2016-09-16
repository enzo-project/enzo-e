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
  array_[0]        = config->initial_sedov_array[0];
  array_[1]        = config->initial_sedov_array[1];
  array_[2]        = config->initial_sedov_array[2];
  radius_relative_ = config->initial_sedov_radius_relative;
  pressure_in_     = config->initial_sedov_pressure_in;
  pressure_out_    = config->initial_sedov_pressure_out;
  density_         = config->initial_sedov_density;
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
 Block * block,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  ASSERT("EnzoInitialSedovArray3",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  ASSERT("EnzoInitialSedovArray3",
	 "Insufficient number of fields",
	 field.field_count() >= 4);

  enzo_float *  d = (enzo_float *) field.values
    (field.field_id("density"));
  
  enzo_float * te = 0;

  te = (enzo_float *) field.values
    (field.field_id("total_energy"));

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  double xmb,ymb,zmb;
  block->data()->lower(&xmb,&ymb,&zmb);

  double xpb,ypb,zpb;
  block->data()->upper(&xpb,&ypb,&zpb);

  double hx,hy,hz;
  field.cell_width(xmb,xpb,&hx,
		   ymb,ypb,&hy,
		   zmb,zpb,&hz);

  // Parameters

  const double sedov_radius = radius_relative_/array_[0];
  const double sedov_radius_2 = sedov_radius*sedov_radius;

  const int in = cello::index_static();

  const double sedov_te_in = 
    pressure_in_  / ((EnzoBlock::Gamma[in] - 1.0) * density_);
  const double sedov_te_out= 
    pressure_out_ / ((EnzoBlock::Gamma[in] - 1.0) * density_);

  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // clear all fields

  for (int iv=0; iv<field.field_count(); iv++) {

    enzo_float * array = (enzo_float *) field.values (iv);

    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int i = INDEX(ix,iy,iz,mx,my);
	  array[i] = 0.0;
	}
      }
    }
  }

  // background 

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {

	int i = INDEX(ix,iy,iz,mx,my);
	d[i]  = density_;
	te[i] = sedov_te_out;

      }
    }
  }

  // array of explosions

  double xmd,ymd,zmd;
  hierarchy->lower(&xmd,&ymd,&zmd);
  double xpd,ypd,zpd;
  hierarchy->upper(&xpd,&ypd,&zpd);

  // bounds of possible explosions intersecting this Block

  double r = sedov_radius;
  
  int kxm = MAX((int)floor((xmb-xmd-r)/(xpd-xmd)*array_[0])-1,0);
  int kym = MAX((int)floor((ymb-ymd-r)/(ypd-ymd)*array_[1])-1,0);
  int kzm = MAX((int)floor((zmb-zmd-r)/(zpd-zmd)*array_[2])-1,0);
  int kxp = MIN( (int)ceil((xpb-xmd+r)/(xpd-xmd)*array_[0])+1,array_[0]);
  int kyp = MIN( (int)ceil((ypb-ymd+r)/(ypd-ymd)*array_[1])+1,array_[1]);
  int kzp = MIN( (int)ceil((zpb-zmd+r)/(zpd-zmd)*array_[2])+1,array_[2]);
  
  TRACE3 ("SEDOV: %d %d %d",kxp,kyp,kzp);

  double hxa = (xpd-xmd) / array_[0];
  double hya = (ypd-ymd) / array_[1];
  double hza = (zpd-zmd) / array_[2];

  // (kx,ky,kz) index bounds of explosions in domain

  for (int kz=kzm; kz<kzp; kz++) {
    double zc = hza*(0.5+kz);
    for (int ky=kym; ky<kyp; ky++) {
      double yc = hya*(0.5+ky);
      for (int kx=kxm; kx<kxp; kx++) {
	double xc = hxa*(0.5+kx);
	TRACE3("xc,yc,zc = %f %f %f",xc,yc,zc);

	// (explosion center xc,yc,zc)

	for (int iz=0; iz<mz; iz++) {
	  double z = zmb + (iz - gz + 0.5)*hz - zc;
	  for (int iy=0; iy<my; iy++) {
	    double y = ymb + (iy - gy + 0.5)*hy - yc;
	    for (int ix=0; ix<mx; ix++) {
	      double x = xmb + (ix - gx + 0.5)*hx - xc;
	      double r2 = x*x + y*y + z*z;

	      int i = INDEX(ix,iy,iz,mx,my);
	  
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

