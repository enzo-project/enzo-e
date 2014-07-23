// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialTurbulence.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:30:49 UTC 2014
/// @brief    Implementation of Enzo 2D Implosion problem initialization

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialTurbulence::EnzoInitialTurbulence 
(int init_cycle, double init_time) throw ()
  : Initial(init_cycle, init_time) 
{ }

//----------------------------------------------------------------------

void EnzoInitialTurbulence::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);
}

//----------------------------------------------------------------------

void EnzoInitialTurbulence::enforce_block 
(
 CommBlock * comm_block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  INCOMPLETE("EnzoInitialTurbulence::enforce_block()");

  ASSERT("EnzoInitialTurbulence",
	 "CommBlock does not exist",
	 comm_block != NULL);

  FieldBlock * field_block = comm_block->block()->field_block();

  enzo_float *  d = (enzo_float *) field_block->values("density");
  enzo_float * ax = (enzo_float *) field_block->values("acceleration_x");
  enzo_float * ay = (enzo_float *) field_block->values("acceleration_y");
  enzo_float * az = (enzo_float *) field_block->values("acceleration_z");
  enzo_float * vx = (enzo_float *) field_block->values("velocity_x");
  enzo_float * vy = (enzo_float *) field_block->values("velocity_y");
  enzo_float * vz = (enzo_float *) field_block->values("velocity_z");
  enzo_float * te = (enzo_float *) field_block->values("total_energy");

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  double xm,ym,zm;
  comm_block->block()->lower(&xm,&ym,&zm);

  double xp,yp,zp;
  comm_block->block()->upper(&xp,&yp,&zp);

  double hx,hy,hz;
  field_block->cell_width(xm,xp,&hx,
			  ym,yp,&hy,
			  zm,zp,&hz);

  int gx,gy,gz;
  field_descr->ghosts(0,&gx,&gy,&gz);
  WARNING("EnzoInitialTurbulence",
	  "Assumes same ghost zone depth for all fields");

  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;
  int ndz = nz + 2*gz;

  for (int iz=gz; iz<nz+gz; iz++) {
    double z = zm + (iz - gz + 0.5)*hz;
    for (int iy=gy; iy<ny+gy; iy++) {
      double y = ym + (iy - gy + 0.5)*hy;
      for (int ix=gx; ix<nx+gx; ix++) {
	double x = xm + (ix - gx + 0.5)*hx;
	int i = ix + ndx*(iy + ndy*iz);
	d[i]  = 0.125;
	vx[i] = 0.0;
	if (vy) vy[i] = 0.0;
	if (vz) vz[i] = 0.0;
	ax[i] = 0.0;
	if (ay) ay[i] = 0.0;
	if (az) az[i] = 0.0;
	te[i] = 0.14 / ((EnzoBlock::Gamma - 1.0) * d[i]);
      }
    }
  }
}
