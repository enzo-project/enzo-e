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

#ifdef CONFIG_USE_CHARM

void EnzoInitialImplosion2::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);
}
#endif

//----------------------------------------------------------------------

void EnzoInitialImplosion2::enforce 
(
 Block * block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  ASSERT("EnzoInitialImplosion2",
	 "Block does not exist",
	 block != NULL);

  FieldBlock * field_block = block->field_block();

  ASSERT("EnzoInitialImplosion2",
	 "Insufficient number of fields",
	 field_descr->field_count() >= 4);

  WARNING("EnzoInitialImplosion2::enforce",
	  "hard-coded field index ordering");

  int index_density         = 0;
  int index_velocity_x      = 1;
  int index_velocity_y      = 2;
  int index_total_energy    = 3;
 
  enzo_float *  d = (enzo_float *) field_block->field_values(index_density);
  enzo_float * vx = (enzo_float *) field_block->field_values(index_velocity_x);
  enzo_float * vy = (enzo_float *) field_block->field_values(index_velocity_y);
  enzo_float * te = (enzo_float *) field_block->field_values(index_total_energy);

  // Block size (excluding ghosts)
  int nx,ny;
  field_block->size(&nx,&ny);

  // Cell widths
  double xm,ym;
  block->lower(&xm,&ym);

  double xp,yp;
  block->upper(&xp,&yp);

  double hx,hy;
  field_block->cell_width(xm,xp,&hx,ym,yp,&hy);

  // Ghost depths
  int gx,gy;
  field_descr->ghosts(index_density,&gx,&gy);

  // Left edges
  double x0, y0;
  block->lower(&x0,&y0);

  // WARNING("EnzoInitialImplosion2",
  // 		  "Assumes same ghost zone depth for all fields");

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;

  for (int iy=gy; iy<ny+gy; iy++) {
    double y = y0 + (iy - gy + 0.5)*hy;
    for (int ix=gy; ix<nx+gx; ix++) {
      double x = x0 + (ix - gx + 0.5)*hx;
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
