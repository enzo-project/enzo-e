// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of Enzo 2D Implosion problem initialization

#include "cello.hpp"

#include "enzo.hpp"

#ifdef CONFIG_USE_GRACKLE

//----------------------------------------------------------------------

EnzoInitialGrackleTest::EnzoInitialGrackleTest 
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time),
    chemistry_(& config->method_grackle_chemistry),
    units_(& config->method_grackle_units)
{ 
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::enforce_block 
(
 CommBlock * comm_block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  ASSERT("EnzoInitialGrackleTest",
	 "CommBlock does not exist",
	 comm_block != NULL);

  FieldBlock * field_block = comm_block->block()->field_block();


  gr_float * density 
    = (gr_float *) field_block->values("density");
  gr_float * internal_energy  
    = (gr_float *) field_block->values("internal_energy");
  gr_float * total_energy  
    = (gr_float *) field_block->values("total_energy");
  gr_float * velocity_x
    = (gr_float *) field_block->values("velocity_x");
  gr_float * velocity_y
    = (gr_float *) field_block->values("velocity_y");
  gr_float * velocity_z
    = (gr_float *) field_block->values("velocity_z");
  gr_float * HI_density 
    = (gr_float *) field_block->values("HI_density");
  gr_float * HII_density 
    = (gr_float *) field_block->values("HII_density");
  gr_float * HM_density 
    = (gr_float *) field_block->values("HM_density");
  gr_float * HeI_density 
    = (gr_float *) field_block->values("HeI_density");
  gr_float * HeII_density 
    = (gr_float *) field_block->values("HeII_density");
  gr_float * HeIII_density 
    = (gr_float *) field_block->values("HeIII_density");
  gr_float * H2I_density 
    = (gr_float *) field_block->values("H2I_density");
  gr_float * H2II_density 
    = (gr_float *) field_block->values("H2II_density");
  gr_float * DI_density 
    = (gr_float *) field_block->values("DI_density");
  gr_float * DII_density 
    = (gr_float *) field_block->values("DII_density");
  gr_float * HDI_density 
    = (gr_float *) field_block->values("HDI_density");
  gr_float * e_density 
    = (gr_float *) field_block->values("e_density");
  gr_float * metal_density 
    = (gr_float *) field_block->values("metal_density");
  gr_float * cooling_time 
    = (gr_float *) field_block->values("cooling_time");
  gr_float * temperature 
    = (gr_float *) field_block->values("temperature");
  gr_float * pressure 
    = (gr_float *) field_block->values("pressure");
  gr_float * gamma 
    = (gr_float *) field_block->values("gamma");

  // Block size (excluding ghosts)
  int nx,ny;
  field_block->size(&nx,&ny);

  // Cell widths
  double xm,ym;
  comm_block->block()->lower(&xm,&ym);

  double xp,yp;
  comm_block->block()->upper(&xp,&yp);

  double hx,hy;
  field_block->cell_width(xm,xp,&hx,ym,yp,&hy);

  // Ghost depths
  int gx,gy;
  field_descr->ghosts(0,&gx,&gy);

  // WARNING("EnzoInitialGrackleTest",
  // 		  "Assumes same ghost zone depth for all fields");

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;

  const double mh     = 1.67262171e-24;
  const double kboltz = 1.3806504e-16;

  gr_float temperature_units =  mh * pow(units_->a_units * 
                                         units_->length_units /
                                         units_->time_units, 2) / kboltz;

  const double density_out = 1.0;
  const double density_in  = 0.125;
  const double pressure_out = 1.0;
  const double pressure_in  = 0.14;
  
  double tiny_number = 1e-20;
  for (int iy=gy; iy<ny+gy; iy++) {
    double y = ym + (iy - gy + 0.5)*hy;
    for (int ix=gx; ix<nx+gx; ix++) {
      double x = xm + (ix - gx + 0.5)*hx;
      int i = INDEX(ix,iy,0,ngx,ngy);

      const bool in = (x + y < 0.3) ;

      density[i] = in ? density_in : density_out;
      internal_energy[i] = in ? (pressure_in) / (0.4 * density[i]) 
	:                       (pressure_out) / (0.4 * density[i]) ;

      HI_density[i]    = density[i] * chemistry_->HydrogenFractionByMass;
      HII_density[i]   = density[i] * tiny_number;
      HM_density[i]    = density[i] * tiny_number;
      HeI_density[i]   = density[i] * (1.0 - chemistry_->HydrogenFractionByMass);
      HeII_density[i]  = density[i] * tiny_number;
      HeIII_density[i] = density[i] * tiny_number;
      H2I_density[i]   = density[i] * tiny_number;
      H2II_density[i]  = density[i] * tiny_number;
      DI_density[i]    = density[i] * 2.0 * 3.4e-5;
      DII_density[i]   = density[i] * tiny_number;
      HDI_density[i]   = density[i] * tiny_number;
      e_density[i]     = density[i] * tiny_number;
      // solar metallicity
      metal_density[i] = density[i] * chemistry_->SolarMetalFractionByMass;

      velocity_x[i] = 0.0;
      velocity_y[i] = 0.0;
      velocity_z[i] = 0.0;

    }
  }
}

#endif /* CONFIG_USE_GRACKLE */
