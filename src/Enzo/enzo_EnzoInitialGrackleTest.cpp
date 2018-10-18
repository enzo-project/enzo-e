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

  grackle_fields_.HI_density      = field.is_field("HI_density") ?
                       (gr_float *) field.values("HI_density")     : NULL;
  grackle_fields_.HII_density     = field.is_field("HII_density") ?
                       (gr_float *) field.values("HII_density")    : NULL;
  grackle_fields_.HM_density      = field.is_field("HM_density") ?
                       (gr_float *) field.values("HM_density")     : NULL;
  grackle_fields_.HeI_density     = field.is_field("HeI_density") ?
                       (gr_float *) field.values("HeI_density")    : NULL;
  grackle_fields_.HeII_density    = field.is_field("HeII_density") ?
                       (gr_float *) field.values("HeII_density")   : NULL;
  grackle_fields_.HeIII_density   = field.is_field("HeIII_density") ?
                       (gr_float *) field.values("HeIII_density")  : NULL;
  grackle_fields_.e_density       = field.is_field("e_density") ?
                       (gr_float *) field.values("e_density")      : NULL;


  grackle_fields_.H2I_density     = field.is_field("H2I_density") ?
                       (gr_float *) field.values("H2I_density") : NULL;
  grackle_fields_.H2II_density    = field.is_field("H2II_density") ?
                       (gr_float *) field.values("H2II_density") : NULL;
  grackle_fields_.DI_density      = field.is_field("DI_density") ?
                       (gr_float *) field.values("DI_density") : NULL;
  grackle_fields_.DII_density     = field.is_field("DII_density") ?
                       (gr_float *) field.values("DII_density") : NULL;
  grackle_fields_.HDI_density     = field.is_field("HDI_Density") ?
                       (gr_float *) field.values("HDI_density") : NULL;
  grackle_fields_.metal_density   = field.is_field("metal_density") ?
                       (gr_float *) field.values("metal_density") : NULL;
  //grackle_fields_.cooling_time  = (gr_float *) field.values("cooling_time");
  //grackle_fields_.temperature   = (gr_float *) field.values("temperature");
  //grackle_fields_.pressure      = (gr_float *) field.values("pressure");
  //grackle_fields_.gamma         = (gr_float *) field.values("gamma");

  gr_float * total_energy  = (gr_float *) field.values("total_energy");
  gr_float * gamma         = (gr_float *) field.values("gamma");
  enzo_float * pressure    = field.is_field("pressure") ?
               (enzo_float*) field.values("pressure") : NULL;
  enzo_float * temperature = field.is_field("temperature") ?
               (enzo_float*) field.values("temperature") : NULL;


  // Block size (excluding ghosts)
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // Cell widths
  double xm,ym,zm;
  block->data()->lower(&xm,&ym,&zm);

  double xp,yp,zp;
  block->data()->upper(&xp,&yp,&zp);

  double hx,hy,hz;
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  // Ghost depths
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  // WARNING("EnzoInitialGrackleTest",
  //      "Assumes same ghost zone depth for all fields");

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  double a_units = 1.0 / (1.0 + enzo_config->physics_cosmology_initial_redshift);

  const double mh = 1.67262171E-24;
  const double kboltz = 1.3806504E-16;

  gr_float temperature_units =  mh * pow(a_units *
                                         enzo_units->velocity(), 2) / kboltz;

  double H_n_slope = log10(enzo_config->initial_grackle_test_maximum_H_number_density /
                           enzo_config->initial_grackle_test_minimum_H_number_density) /
                           double(nx);
  double metallicity_slope = log10(enzo_config->initial_grackle_test_maximum_metallicity/
                                   enzo_config->initial_grackle_test_minimum_metallicity)/
                                   double(ny);
  double temperature_slope = log10(enzo_config->initial_grackle_test_maximum_temperature/
                                   enzo_config->initial_grackle_test_minimum_temperature)/
                                   double(nz);


  double tiny_number = 1e-20;

  for (int iz=gz; iz<nz+gz; iz++){ // Temperature
    for (int iy=gy; iy<ny+gy; iy++) { // Metallicity
      for (int ix=gx; ix<nx+gx; ix++) { // H Number Density
        int i = INDEX(ix,iy,iz,ngx,ngy);

// AE NEED TO FIX
        grackle_fields_.density[i] = cello::mass_hydrogen *
                     pow(10.0, ((H_n_slope * (ix-gx)) + log10(enzo_config->initial_grackle_test_minimum_H_number_density)))/
                     grackle_data->HydrogenFractionByMass / enzo_units->density();
        // solar metallicity
        grackle_fields_.metal_density[i] = pow(10.0,((metallicity_slope * (iy-gy)) +
                                          log10(enzo_config->initial_grackle_test_minimum_metallicity))) *
                                          grackle_data->SolarMetalFractionByMass * grackle_fields_.density[i];
        grackle_fields_.x_velocity[i] = 0.0;
        grackle_fields_.y_velocity[i] = 0.0;
        grackle_fields_.z_velocity[i] = 0.0;

        if (grackle_data->primordial_chemistry > 0){
            grackle_fields_.HI_density[i]    = grackle_fields_.density[i] * grackle_data->HydrogenFractionByMass;
            grackle_fields_.HII_density[i]   = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.HeI_density[i]   = grackle_fields_.density[i] * (1.0 - grackle_data->HydrogenFractionByMass);
            grackle_fields_.HeII_density[i]  = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.HeIII_density[i] = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.e_density[i]     = grackle_fields_.density[i] * tiny_number;
        }
        if (grackle_data->primordial_chemistry > 1){
          grackle_fields_.HM_density[i]    = grackle_fields_.density[i] * tiny_number;
          grackle_fields_.H2I_density[i]   = grackle_fields_.density[i] * tiny_number;
          grackle_fields_.H2II_density[i]  = grackle_fields_.density[i] * tiny_number;
        }
        if (grackle_data->primordial_chemistry > 2){
            grackle_fields_.DI_density[i]    = grackle_fields_.density[i] * grackle_data->DeuteriumToHydrogenRatio;
            grackle_fields_.DII_density[i]   = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.HDI_density[i]   = grackle_fields_.density[i] * tiny_number;
        }

      }
    }
  }

  /* Set internal energy and temperature */
  enzo_float mu = enzo_config->ppm_mol_weight;

  for (int iz=gz; iz<nz+gz; iz++){ // Temperature
    for (int iy=gy; iy<ny+gy; iy++) { // Metallicity
      for (int ix=gx; ix<nx+gx; ix++) { // H Number Density
        int i = INDEX(ix,iy,iz,ngx,ngy);

        /* calculate mu if following species */
        if (grackle_data->primordial_chemistry > 0){

          mu = grackle_fields_.density[i] + grackle_fields_.HI_density[i] +
                      grackle_fields_.HeII_density[i] +
                      (grackle_fields_.HeI_density[i] +
                       grackle_fields_.HeII_density[i] +
                       grackle_fields_.HeIII_density[i])*0.25;

          if (grackle_data->primordial_chemistry > 1){

            mu += (grackle_fields_.HM_density[i] + grackle_fields_.H2I_density[i] +
                   grackle_fields_.H2II_density[i])*0.5;
          }
          if (grackle_data->primordial_chemistry > 2){
            mu += (grackle_fields_.DI_density[i] + grackle_fields_.DII_density[i])*0.5 +
                   grackle_fields_.HDI_density[i]/3.0;
          }

          mu = grackle_fields_.density[i] / mu;
        } // end primordial_chemistry > 0

        grackle_fields_.internal_energy[i] = pow(10.0, ((temperature_slope * (iz-gz)) +
                                      log10(enzo_config->initial_grackle_test_minimum_temperature)))/
                             mu / temperature_units / (enzo_config->field_gamma - 1.0);
        total_energy[i]    = grackle_fields_.internal_energy[i];
        gamma[i]           = enzo_config->field_gamma;
      }
    }
  }

  // Finally compute temperature and pressure if fields are tracked
  // for output
  const int in = cello::index_static();
  int comoving_coordinates = enzo_config->physics_cosmology;

  if (pressure){
    EnzoComputePressure compute_pressure (EnzoBlock::Gamma[in],
                                          comoving_coordinates);
    compute_pressure.compute(enzo_block);
  }

  if (temperature){
    EnzoComputeTemperature compute_temperature
      (enzo_config->ppm_density_floor,
       enzo_config->ppm_temperature_floor,
       enzo_config->ppm_mol_weight,
       comoving_coordinates);

    compute_temperature.compute(enzo_block);
  }


#endif /* CONFIG_USE_GRACKLE */
}
