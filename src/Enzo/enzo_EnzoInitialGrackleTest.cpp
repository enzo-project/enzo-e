// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     Tue May 7 2019
/// @brief    [\ref Enzo] Grackle chemistry/cooling library initial conditions

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialGrackleTest::EnzoInitialGrackleTest
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time)
{
  return;
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

#ifdef CONFIG_USE_GRACKLE
  static bool warn[CONFIG_NODE_SIZE] = {false};
  const int in = cello::index_static();
  if (! warn[in]) {
    WARNING("EnzoInitialGrackleTest::pup()","Skipping units_");
    //  const code_units      * units_;
    WARNING("EnzoInitialGrackleTest::pup()", "Skipping chemistry_");
    //  const chemistry_data  * chemistry_;
    warn[in] = true;
  }
#endif
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::enforce_block
(
 Block * block, const Hierarchy  * hierarchy) throw()
{
#ifndef CONFIG_USE_GRACKLE

  ERROR("EnzoInitialGrackleTest::compute()",
  "Trying to use Initialization 'grackle_test' with "
  "Grackle configuration turned off!");

#else /* CONFIG_USE_GRACKLE */

  ASSERT("EnzoInitialGrackleTest",
   "Block does not exist",
   block != NULL);

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits  * enzo_units = enzo::units();

  Field field = block->data()->field();

  grackle_field_data grackle_fields_;

  EnzoMethodGrackle::setup_grackle_fields(enzo_block, & grackle_fields_);

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

  double H_n_slope = log10(enzo_config->initial_grackle_test_maximum_H_number_density /
                           enzo_config->initial_grackle_test_minimum_H_number_density) /
                           double(nx);

  double temperature_slope = log10(enzo_config->initial_grackle_test_maximum_temperature/
                                   enzo_config->initial_grackle_test_minimum_temperature)/
                                   double(ny);

  double metallicity_slope = log10(enzo_config->initial_grackle_test_maximum_metallicity/
                                   enzo_config->initial_grackle_test_minimum_metallicity)/
                                   double(nz);

  double tiny_number = 1e-10;

  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

  for (int iz=0; iz<nz+gz; iz++){ // Metallicity
    for (int iy=0; iy<ny+gy; iy++) { // Temperature
      for (int ix=0; ix<nx+gx; ix++) { // H Number Density
        int i = INDEX(ix,iy,iz,ngx,ngy);

        grackle_fields_.density[i] = cello::mass_hydrogen *
                     pow(10.0, ((H_n_slope * (ix-gx)) + log10(enzo_config->initial_grackle_test_minimum_H_number_density)))/
                     grackle_chemistry->HydrogenFractionByMass / enzo_units->density();

        // solar metallicity
        grackle_fields_.metal_density[i] = pow(10.0,((metallicity_slope * (iz-gz)) +
                                          log10(enzo_config->initial_grackle_test_minimum_metallicity))) *
                                          grackle_chemistry->SolarMetalFractionByMass * grackle_fields_.density[i];
        grackle_fields_.x_velocity[i] = 0.0;
        grackle_fields_.y_velocity[i] = 0.0;
        grackle_fields_.z_velocity[i] = 0.0;

        if (grackle_chemistry->primordial_chemistry > 0){
            grackle_fields_.HI_density[i]    = grackle_fields_.density[i] * grackle_chemistry->HydrogenFractionByMass;
            grackle_fields_.HII_density[i]   = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.HeI_density[i]   = grackle_fields_.density[i] * (1.0 - grackle_chemistry->HydrogenFractionByMass);
            grackle_fields_.HeII_density[i]  = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.HeIII_density[i] = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.e_density[i]     = grackle_fields_.density[i] * tiny_number;
        }
        if (grackle_chemistry->primordial_chemistry > 1){
          grackle_fields_.HM_density[i]    = grackle_fields_.density[i] * tiny_number;
          grackle_fields_.H2I_density[i]   = grackle_fields_.density[i] * tiny_number;
          grackle_fields_.H2II_density[i]  = grackle_fields_.density[i] * tiny_number;
        }
        if (grackle_chemistry->primordial_chemistry > 2){
            grackle_fields_.DI_density[i]    = grackle_fields_.density[i] * grackle_chemistry->DeuteriumToHydrogenRatio;
            grackle_fields_.DII_density[i]   = grackle_fields_.density[i] * tiny_number;
            grackle_fields_.HDI_density[i]   = grackle_fields_.density[i] * tiny_number;
        }

      }
    }
  }

  /* Set internal energy and temperature */
  enzo_float mu = enzo_config->ppm_mol_weight;

  for (int iz=0; iz<nz+gz; iz++){ // Metallicity
    for (int iy=0; iy<ny+gy; iy++) { // Temperature
      for (int ix=0; ix<nx+gx; ix++) { // H Number Density
        int i = INDEX(ix,iy,iz,ngx,ngy);

        /* calculate mu if following species */
        if (grackle_chemistry->primordial_chemistry > 0){

          mu = grackle_fields_.e_density[i] + grackle_fields_.HI_density[i] +
                      grackle_fields_.HII_density[i] +
                      (grackle_fields_.HeI_density[i] +
                       grackle_fields_.HeII_density[i] +
                       grackle_fields_.HeIII_density[i])*0.25;

          if (grackle_chemistry->primordial_chemistry > 1){

            mu += grackle_fields_.HM_density[i] + (grackle_fields_.H2I_density[i] +
                   grackle_fields_.H2II_density[i])*0.5;
          }
          if (grackle_chemistry->primordial_chemistry > 2){
            mu += (grackle_fields_.DI_density[i] + grackle_fields_.DII_density[i])*0.5 +
                   grackle_fields_.HDI_density[i]/3.0;
          }

          mu = grackle_fields_.density[i] / mu;
        } // end primordial_chemistry > 0

        grackle_fields_.internal_energy[i] = pow(10.0, ((temperature_slope * (iy-gy)) +
                                      log10(enzo_config->initial_grackle_test_minimum_temperature)))/
                             mu / enzo_units->temperature() / (enzo_config->field_gamma - 1.0);
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

    // Note: using compute_ method to avoid re-generating grackle_fields
    //       struct. Otherwise, one could just do:
    //             compute_pressure.compute(enzo_block, pressure);
    //       OR
    //             compute_pressure.compute(enzo_block);
    //
    //       The former provides ability to compute pressure into an
    //       array that is non-static (i.e. not a field that persists).
    //
    compute_pressure.compute_(enzo_block,
                              pressure,
                              NULL,
                              &grackle_fields_);
  }

  if (temperature){
    EnzoComputeTemperature compute_temperature
      (enzo_config->ppm_density_floor,
       enzo_config->ppm_temperature_floor,
       enzo_config->ppm_mol_weight,
       comoving_coordinates);

    compute_temperature.compute_(enzo_block,
                                 temperature,
                                 false, // do not re-compute pressure field
                                 NULL, &grackle_fields_
                                 );
  }

  EnzoMethodGrackle::delete_grackle_fields(&grackle_fields_);

  return;
#endif /* CONFIG_USE_GRACKLE */
}
