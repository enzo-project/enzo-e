// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     Tue May 7 2019
/// @brief    [\ref Enzo] Grackle chemistry/cooling library initial conditions

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/tests/tests.hpp"

//----------------------------------------------------------------------

EnzoInitialGrackleTest::EnzoInitialGrackleTest
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time),
    min_max_H_number_density_{},
    min_max_metallicity_{},
    min_max_temperature_{}
{
  min_max_H_number_density_ = {config->initial_grackle_test_minimum_H_number_density,
                               config->initial_grackle_test_maximum_H_number_density};
  min_max_metallicity_ = {config->initial_grackle_test_minimum_metallicity,
                          config->initial_grackle_test_maximum_metallicity};
  min_max_temperature_ = {config->initial_grackle_test_minimum_temperature,
                          config->initial_grackle_test_maximum_temperature};
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | min_max_H_number_density_;
  p | min_max_metallicity_;
  p | min_max_temperature_;
}

//----------------------------------------------------------------------

void EnzoInitialGrackleTest::enforce_block
(
 Block * block, const Hierarchy  * hierarchy) throw()
{

  ASSERT("EnzoInitialGrackleTest",
   "Block does not exist",
   block != NULL);

  EnzoUnits  * enzo_units = enzo::units();

  Field field = block->data()->field();

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

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;

  double H_n_slope = log10(min_max_H_number_density_[1]/
                           min_max_H_number_density_[0]) / double(nx);

  double temperature_slope = log10(min_max_temperature_[1]/
                                   min_max_temperature_[0]) / double(ny);

  double metallicity_slope = log10(min_max_metallicity_[1]/
                                   min_max_metallicity_[0]) / double(nz);

  double tiny_number = 1e-10;

  // load some information about how grackle is configured
  const GrackleChemistryData * grackle_chem = enzo::grackle_chemistry();
  ASSERT("EnzoInitialGrackleTest",
         "the simulation must be configured to use grackle",
         grackle_chem != nullptr);
  const int primordial_chemistry
    = grackle_chem->get<int>("primordial_chemistry");

  // load some constants from the Grackle configuration
  const double HydrogenFractionByMass
    = grackle_chem->get<double>("HydrogenFractionByMass");
  const double SolarMetalFractionByMass
    = grackle_chem->get<double>("SolarMetalFractionByMass");
  const double DeuteriumToHydrogenRatio
    = grackle_chem->get<double>("DeuteriumToHydrogenRatio");

  // load in fields required in all grackle-configurations
  CelloView<enzo_float,3> density    = field.view<enzo_float>("density");
  CelloView<enzo_float,3> velocity_x = field.view<enzo_float>("velocity_x");
  CelloView<enzo_float,3> velocity_y = field.view<enzo_float>("velocity_y");
  CelloView<enzo_float,3> velocity_z = field.view<enzo_float>("velocity_z");
  CelloView<enzo_float,3> eint = field.view<enzo_float>("internal_energy");
  CelloView<enzo_float,3> etot = field.view<enzo_float>("total_energy");

  // at the time of writing the following check, it looks like grackle
  // will soon unlock extra functionality with metal_cooling > 1
  ASSERT("EnzoInitialGrackleTest",
         "the initializer currently assumes Grackle has metal_cooling=1",
         grackle_chem->get<int>("metal_cooling") == 1);
  CelloView<enzo_float,3> metal_density = field.view<enzo_float>("metal_density");

  // now we load in species fields
  CelloView<enzo_float,3> HI_density = (primordial_chemistry > 0) ?
    field.view<enzo_float>("HI_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> HII_density = (primordial_chemistry > 0) ?
    field.view<enzo_float>("HII_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> HeI_density = (primordial_chemistry > 0) ?
    field.view<enzo_float>("HeI_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> HeII_density = (primordial_chemistry > 0) ?
    field.view<enzo_float>("HeII_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> HeIII_density = (primordial_chemistry > 0) ?
    field.view<enzo_float>("HeIII_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> e_density = (primordial_chemistry > 0) ?
    field.view<enzo_float>("e_density") : CelloView<enzo_float,3>();

  CelloView<enzo_float,3> HM_density = (primordial_chemistry > 1) ?
    field.view<enzo_float>("HM_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> H2I_density = (primordial_chemistry > 1) ?
    field.view<enzo_float>("H2I_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> H2II_density = (primordial_chemistry > 1) ?
    field.view<enzo_float>("H2II_density") : CelloView<enzo_float,3>();

  CelloView<enzo_float,3> DI_density = (primordial_chemistry > 2) ?
    field.view<enzo_float>("DI_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> DII_density = (primordial_chemistry > 2) ?
    field.view<enzo_float>("DII_density") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> HDI_density = (primordial_chemistry > 2) ?
    field.view<enzo_float>("HDI_density") : CelloView<enzo_float,3>(); 

  for (int iz=0; iz<nz+gz; iz++){ // Metallicity
    for (int iy=0; iy<ny+gy; iy++) { // Temperature
      for (int ix=0; ix<nx+gx; ix++) { // H Number Density

        const enzo_float cur_density = 
          (enzo_constants::mass_hydrogen *
           pow(10.0, ((H_n_slope * (ix-gx)) + log10(min_max_H_number_density_[0])))/
           HydrogenFractionByMass / enzo_units->density());

        density(iz,iy,ix) = cur_density;

        // solar metallicity
        enzo_float solar_metal_dens = SolarMetalFractionByMass * cur_density;
        metal_density(iz,iy,ix) = pow(10.0,((metallicity_slope * (iz-gz)) +
                                            log10(min_max_metallicity_[0]))) *
                                  solar_metal_dens;
        velocity_x(iz,iy,ix) = 0.0;
        velocity_y(iz,iy,ix) = 0.0;
        velocity_z(iz,iy,ix) = 0.0;

        if (primordial_chemistry > 0){
          HI_density(iz,iy,ix)    = cur_density * HydrogenFractionByMass;
          HII_density(iz,iy,ix)   = cur_density * tiny_number;
          HeI_density(iz,iy,ix)   = cur_density * (1.0 - HydrogenFractionByMass);
          HeII_density(iz,iy,ix)  = cur_density * tiny_number;
          HeIII_density(iz,iy,ix) = cur_density * tiny_number;
          e_density(iz,iy,ix)     = cur_density * tiny_number;
        }
        if (primordial_chemistry > 1){
          HM_density(iz,iy,ix)    = cur_density * tiny_number;
          H2I_density(iz,iy,ix)   = cur_density * tiny_number;
          H2II_density(iz,iy,ix)  = cur_density * tiny_number;
        }
        if (primordial_chemistry > 2){
          DI_density(iz,iy,ix)    = cur_density * DeuteriumToHydrogenRatio;
          DII_density(iz,iy,ix)   = cur_density * tiny_number;
          HDI_density(iz,iy,ix)   = cur_density * tiny_number;
        }

      }
    }
  }

  /* Set internal energy from temperature requirements */
  enzo_float nominal_mu = enzo::fluid_props()->mol_weight();
  const enzo_float nominal_gamma = enzo::fluid_props()->gamma();

  for (int iz=0; iz<nz+gz; iz++){ // Metallicity
    for (int iy=0; iy<ny+gy; iy++) { // Temperature
      for (int ix=0; ix<nx+gx; ix++) { // H Number Density
        int i = INDEX(ix,iy,iz,ngx,ngy);

        enzo_float mu;
        if (primordial_chemistry == 0) {
          mu = nominal_mu;
        } else {
          enzo_float tmp = 
            ( (e_density(iz,iy,ix) + HI_density(iz,iy,ix) +
               HII_density(iz,iy,ix)) +
              0.25 * (HeI_density(iz,iy,ix) + HeII_density(iz,iy,ix) +
                      HeIII_density(iz,iy,ix) ) );

          if (primordial_chemistry > 1){
            tmp += HM_density(iz,iy,ix) + 0.5 * (H2I_density(iz,iy,ix) +
                                                 H2II_density(iz,iy,ix));
          }
          if (primordial_chemistry > 2){
            tmp += (0.5 * (DI_density(iz,iy,ix) + DII_density(iz,iy,ix))
                    + HDI_density(iz,iy,ix)/3.0);
          }

          mu = density(iz,iy,ix) / tmp;
        } // end primordial_chemistry > 0

        enzo_float specific_thermal_energy =
          (pow(10.0, ((temperature_slope * (iy-gy)) +
                      log10(min_max_temperature_[0]))) /
           mu / enzo_units->kelvin_per_energy_units() / (nominal_gamma - 1.0));
        eint(iz,iy,ix) = specific_thermal_energy;
        etot(iz,iy,ix) = specific_thermal_energy;
      }
    }
  }

  // In the original version of this method, we computed the temperature and
  // pressure fields here, but there wasn't any point to doing that.
  // -> I think it sets a somewhat poor example to do that without any reason,
  //    so I have deleted that logic

  block->initial_done();

  return;
}
