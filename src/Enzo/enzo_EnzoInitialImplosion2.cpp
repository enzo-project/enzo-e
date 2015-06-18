// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialImplosion2.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of Enzo 2D Implosion problem initialization

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialImplosion2::EnzoInitialImplosion2 
(int init_cycle, double init_time) throw ()
  : Initial(init_cycle, init_time) 
{ }

//----------------------------------------------------------------------

void EnzoInitialImplosion2::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);
}

//----------------------------------------------------------------------

void EnzoInitialImplosion2::enforce_block 
(
 Block * block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  ASSERT("EnzoInitialImplosion2",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  ASSERT("EnzoInitialImplosion2::enforce_block",
	 "Insufficient number of fields",
	 field.field_count() >= 4);

  WARNING("EnzoInitialImplosion2::enforce_block",
	  "hard-coded field index ordering");

  enzo_float *  d = (enzo_float *) field.values("density");
  enzo_float * vx = (enzo_float *) field.values("velocity_x");
  enzo_float * vy = (enzo_float *) field.values("velocity_y");
  enzo_float * te = (enzo_float *) field.values("total_energy");

  // Block size (excluding ghosts)
  int nx,ny;
  field.size(&nx,&ny);

  // Cell widths
  double xm,ym;
  block->data()->lower(&xm,&ym);

  double xp,yp;
  block->data()->upper(&xp,&yp);

  double hx,hy;
  field.cell_width(xm,xp,&hx,ym,yp,&hy);

  // Ghost depths
  int gx,gy;
  field.ghost_depth(0,&gx,&gy);

  // WARNING("EnzoInitialImplosion2",
  // 		  "Assumes same ghost zone depth for all fields");

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;

  for (int iy=gy; iy<ny+gy; iy++) {
    double y = ym + (iy - gy + 0.5)*hy;
    for (int ix=gx; ix<nx+gx; ix++) {
      double x = xm + (ix - gx + 0.5)*hx;
      int i = INDEX(ix,iy,0,ngx,ngy);
      if (x + y < 0.1517) {
	d[i]  = 0.125;
	vx[i] = 0.0;
	vy[i] = 0.0;
	te[i] = 0.14 / ((EnzoBlock::Gamma - 1.0) * d[i]);
      } else {
	d[i]  = 1.0;
	vx[i] = 0.0;
	vy[i] = 0.0;
	te[i] = 1.0 / ((EnzoBlock::Gamma - 1.0) * d[i]);
      }
    }
  }

}
