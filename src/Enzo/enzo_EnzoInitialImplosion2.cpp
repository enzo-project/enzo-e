// $Id: enzo_EnzoInitialImplosion2.cpp 1877 2010-11-30 01:20:27Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialImplosion2.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of Enzo 2D Implosion problem initialization

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialImplosion2::EnzoInitialImplosion2
(
 Error     * error,
 Monitor   * monitor,
 EnzoDescr * enzo_descr
 ) throw ()
  : Initial(error, monitor),
    enzo_descr_(enzo_descr)
{}

//----------------------------------------------------------------------

void EnzoInitialImplosion2::initialize (DataDescr * data_descr) throw()
{
  // Intentionally blank
}

//----------------------------------------------------------------------

void EnzoInitialImplosion2::finalize (DataDescr * data_descr) throw()
{
  // Intentionally blank
}

//----------------------------------------------------------------------

void EnzoInitialImplosion2::initialize_block (DataBlock * data_block) throw()

{
  FieldBlock * field_block = data_block->field_block();
  FieldDescr * field_descr = data_block->field_block()->field_descr();

  ASSERT("EnzoInitialImplosion2",
	 "Insufficient number of fields",
	 field_descr->field_count() >= 4);

  int index_density         = 0;
  int index_velocity_x      = 1;
  int index_velocity_y      = 2;
  int index_total_energy    = 3;
 
  Scalar *  d = (Scalar * ) field_block->field_values(index_density);
  Scalar * vx = (Scalar * ) field_block->field_values(index_velocity_x);
  Scalar * vy = (Scalar * ) field_block->field_values(index_velocity_y);
  Scalar * te = (Scalar * ) field_block->field_values(index_total_energy);

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  double hx,hy,hz;
  field_block->cell_width(&hx,&hy,&hz);

 int gx,gy,gz;
  field_descr->ghosts(index_density,&gx,&gy,&gz);

  // WARNING_MESSAGE("EnzoInitialImplosion2",
  // 		  "Assumes same ghost zone depth for all fields");

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;

  for (int iy=gy; iy<ny+gy; iy++) {
    double y = (iy - gy + 0.5)*hy;
    for (int ix=gy; ix<nx+gx; ix++) {
      double x = (ix - gx + 0.5)*hx;
      int i = INDEX(ix,iy,0,ngx,ngy);
      if (x + y < 0.1517) {
	d[i]  = 0.125;
	vx[i] = 0;
	vy[i] = 0;
	te[i] = 0.14 / ((enzo_descr_->Gamma - 1.0) * d[i]);
      } else {
	d[i]  = 1.0;
	vx[i] = 0;
	vy[i] = 0;
	te[i] = 1.0 / ((enzo_descr_->Gamma - 1.0) * d[i]);
      }
    }
  }

}

//----------------------------------------------------------------------

void EnzoInitialImplosion2::finalize_block (DataBlock * data_block) throw()
{
  // Intentionally blank
}

//----------------------------------------------------------------------

std::string EnzoInitialImplosion2::method_name() const throw()
{
  return "implosion_2d";
}

//----------------------------------------------------------------------
