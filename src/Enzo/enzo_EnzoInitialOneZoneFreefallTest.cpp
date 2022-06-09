// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     Tue May 7 2019
/// @brief    [\ref Enzo] Grackle chemistry/cooling library initial conditions

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialOneZoneFreefallTest::EnzoInitialOneZoneFreefallTest
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time)
{
  return;
}

//----------------------------------------------------------------------

void EnzoInitialOneZoneFreefallTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

#ifdef CONFIG_USE_GRACKLE
  static bool warn[CONFIG_NODE_SIZE] = {false};
  const int in = cello::index_static();
  if (! warn[in]) {
    WARNING("EnzoInitialOneZoneFreefallTest::pup()","Skipping units_");
    //  const code_units      * units_;
    WARNING("EnzoInitialOneZoneFreefallTest::pup()", "Skipping chemistry_");
    //  const chemistry_data  * chemistry_;
    warn[in] = true;
  }
#endif
}

//----------------------------------------------------------------------

void EnzoInitialOneZoneFreefallTest::enforce_block
(
 Block * block, const Hierarchy  * hierarchy) throw()
{
#ifndef CONFIG_USE_GRACKLE

  ERROR("EnzoInitialOneZoneFreefallTest::compute()",
  "Trying to use Initialization 'one_zone_freefall_test' with "
  "Grackle configuration turned off!");

#else /* CONFIG_USE_GRACKLE */

  ASSERT("EnzoInitialOneZoneFreefallTest",
   "Block does not exist",
   block != NULL);

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits  * enzo_units = enzo::units();

  Field field = block->data()->field();

  grackle_field_data grackle_fields_;

  const EnzoMethodGrackle * grackle_method = enzo::grackle_method();

  grackle_method->setup_grackle_fields(enzo_block, & grackle_fields_);

  gr_float * total_energy  = (gr_float *) field.values("total_energy");
  gr_float * metal_density = (gr_float *) field.values("metal_density");

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

  double initial_density = enzo_config->initial_one_zone_freefall_test_density;
  double minimum_energy = enzo_config->initial_one_zone_freefall_test_minimum_energy;
  double maximum_energy = enzo_config->initial_one_zone_freefall_test_maximum_energy;

  double minimum_metallicity = enzo_config->initial_one_zone_freefall_test_minimum_metallicity;
  double maximum_metallicity = enzo_config->initial_one_zone_freefall_test_maximum_metallicity;

  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

  // Set up a 1d grid that varies in energy or a 2d grid that also varies in metallicity
  double tiny_number = 1e-10;

  for (int iz=0; iz<nz+gz; iz++){ // nothing
    for (int iy=0; iy<ny+gy; iy++) { // metallicity
      for (int ix=0; ix<nx+gx; ix++) { // energy
        int i = INDEX(ix,iy,iz,ngx,ngy);

        grackle_fields_.density[i] = initial_density;

        // set energy
        total_energy[i] = minimum_energy * pow( maximum_energy/minimum_energy, ix/nx );

        grackle_fields_.internal_energy[i] = total_energy[i];

        // set metallicity
        grackle_fields_.metal_density[i] = minimum_metallicity*cello::metallicity_solar*initial_density *
                                             pow( maximum_metallicity/minimum_metallicity, iy/ny );

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

        // set velocities
        grackle_fields_.x_velocity[i] = 0.0;
        grackle_fields_.y_velocity[i] = 0.0;
        grackle_fields_.z_velocity[i] = 0.0;

      }
    }
  }

  grackle_method->delete_grackle_fields(&grackle_fields_);

  block->initial_done();

  return;
#endif /* CONFIG_USE_GRACKLE */
}
