// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    [\ref Enzo] Grackle chemistry/cooling library initial conditions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialGrackleTest::EnzoInitialGrackleTest
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time)
//  #ifdef CONFIG_USE_GRACKLE
//    , units_()
//  #endif
{
//#ifdef CONFIG_USE_GRACKLE
    //chemistry_ = config->method_grackle_chemistry;
    //units_     = config->method_grackle_units;
//#endif
  return;
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

#ifdef CONFIG_USE_GRACKLE
  //p | units_;
  /*
  WARNING("EnzoInitialGrackleTest::pup()","Skipping units_");
    code_units        units_;
  WARNING("EnzoInitialGrackleTest::pup()", "Skipping chemistry_");
    chemistry_data  * chemistry_;
  */
#endif
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::enforce_block
(
 Block * block,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 const Hierarchy  * hierarchy
 ) throw()

{
#ifndef CONFIG_USE_GRACKLE

  ERROR("EnzoInitialGrackleTest::compute()",
  "Trying to use Initialization 'grackle_test' with "
  "Grackle configuration turned off!");

#else /* CONFIG_USE_GRACKLE */

  ASSERT("EnzoInitialGrackleTest",
   "Block does not exist",
   block != NULL);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  const EnzoConfig * enzo_config = static_cast<const EnzoConfig*>
       (enzo_block->simulation()->config());
  EnzoUnits * enzo_units = (EnzoUnits *) enzo_block->simulation()->problem()->units();
  // units_ = enzo_config->method_grackle_units;

  Field field = block->data()->field();

  grackle_field_data grackle_fields_;

  grackle_fields_.density         = (gr_float *) field.values("density");
  grackle_fields_.internal_energy = (gr_float *) field.values("internal_energy");
  grackle_fields_.x_velocity      = (gr_float *) field.values("velocity_x");
  grackle_fields_.y_velocity      = (gr_float *) field.values("velocity_y");
  grackle_fields_.z_velocity      = (gr_float *) field.values("velocity_z");
  grackle_fields_.HI_density      = (gr_float *) field.values("HI_density");
  grackle_fields_.HII_density     = (gr_float *) field.values("HII_density");
  grackle_fields_.HM_density      = (gr_float *) field.values("HM_density");
  grackle_fields_.HeI_density     = (gr_float *) field.values("HeI_density");
  grackle_fields_.HeII_density    = (gr_float *) field.values("HeII_density");
  grackle_fields_.HeIII_density   = (gr_float *) field.values("HeIII_density");
  grackle_fields_.H2I_density     = (gr_float *) field.values("H2I_density");
  grackle_fields_.H2II_density    = (gr_float *) field.values("H2II_density");
  grackle_fields_.DI_density      = (gr_float *) field.values("DI_density");
  grackle_fields_.DII_density     = (gr_float *) field.values("DII_density");
  grackle_fields_.HDI_density     = (gr_float *) field.values("HDI_density");
  grackle_fields_.e_density       = (gr_float *) field.values("e_density");
  grackle_fields_.metal_density   = (gr_float *) field.values("metal_density");
  //grackle_fields_.cooling_time  = (gr_float *) field.values("cooling_time");
  //grackle_fields_.temperature   = (gr_float *) field.values("temperature");
  //grackle_fields_.pressure      = (gr_float *) field.values("pressure");
  //grackle_fields_.gamma         = (gr_float *) field.values("gamma");

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

  // WARNING("EnzoInitialGrackleTest",
  //      "Assumes same ghost zone depth for all fields");

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;

  const double mh     = 1.67262171e-24;
  const double kboltz = 1.3806504e-16;
  double a_units = 1.0 / (1.0 + enzo_config->physics_cosmology_initial_redshift);

  gr_float temperature_units =  mh * pow(a_units *
                                         enzo_units->length()/
                                         enzo_units->time(), 2) / kboltz;

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

      grackle_fields_.density[i] = in ? density_in : density_out;
      grackle_fields_.internal_energy[i] = in ? (pressure_in)  / (0.4 * grackle_fields_.density[i])
                                              : (pressure_out) / (0.4 * grackle_fields_.density[i]);

      grackle_fields_.HI_density[i]    = grackle_fields_.density[i] * grackle_data->HydrogenFractionByMass;
      grackle_fields_.HII_density[i]   = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.HM_density[i]    = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.HeI_density[i]   = grackle_fields_.density[i] * (1.0 - grackle_data->HydrogenFractionByMass);
      grackle_fields_.HeII_density[i]  = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.HeIII_density[i] = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.H2I_density[i]   = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.H2II_density[i]  = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.DI_density[i]    = grackle_fields_.density[i] * grackle_data->DeuteriumToHydrogenRatio;
      grackle_fields_.DII_density[i]   = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.HDI_density[i]   = grackle_fields_.density[i] * tiny_number;
      grackle_fields_.e_density[i]     = grackle_fields_.density[i] * tiny_number;
      // solar metallicity
      grackle_fields_.metal_density[i] = grackle_fields_.density[i] * grackle_data->SolarMetalFractionByMass;

      grackle_fields_.x_velocity[i] = 0.0;
      grackle_fields_.y_velocity[i] = 0.0;
      grackle_fields_.z_velocity[i] = 0.0;

    }
  }
#endif /* CONFIG_USE_GRACKLE */
}
