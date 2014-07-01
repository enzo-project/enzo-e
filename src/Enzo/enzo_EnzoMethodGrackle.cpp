// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodGrackle class
///
/// Parameters
///
/// code_units units;
///  units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
///  units.density_units = 1.67e-24;
///  units.length_units = 1.0;
///  units.time_units = 1.0e12;
///  units.velocity_units = units.length_units / units.time_units;
///  units.a_units = 1.0; // units for the expansion factor
///
///  chemistry_data chemistry = set_default_chemistry_parameters();
///  chemistry_.use_grackle = 1;            // chemistry on
///  chemistry_.with_radiative_cooling = 1; // cooling on
///  chemistry_.primordial_chemistry = 3;   // molecular network with H, He, D
///  chemistry_.metal_cooling = 1;          // metal cooling on
///  chemistry_.UVbackground = 1;           // UV background on
///  chemistry_.grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5";
///
///  Fields
///
///  density
///  energy
///  x_velocity
///  y_velocity
///  z_velocity
///
/// for primordial_chemistry >= 1
///  HI_density
///  HII_density
///  HeI_density
///  HeII_density
///  HeIII_density
///  e_density
///
///  for primordial_chemistry >= 2
///  HM_density    = new gr_float[field_size];
///  H2I_density   = new gr_float[field_size];
///  H2II_density  = new gr_float[field_size];
///
///  for primordial_chemistry >= 3
///  DI_density    = new gr_float[field_size];
///  DII_density   = new gr_float[field_size];
///  HDI_density   = new gr_float[field_size];
///
/// for metal_cooling = 1
///  metal_density = new gr_float[field_size];

/// Set expansion factor to 1 for non-cosmological simulation.
///  gr_float initial_redshift = 100.;
///  gr_float a_value = 1. / (1. + initial_redshift);

  // Finally, initialize the chemistry object.
///  if (initialize_chemistry_data(chemistry, units, a_value) == 0) {
///    fprintf(stderr, "Error in initialize_chemistry_data.\n");
///    return 0;
//  }

  // Allocate field arrays.
//  gr_float tiny_number = 1.e-20;
//
//  // Set grid dimension and size.
//  // grid_start and grid_end are used to ignore ghost zones.
//  gr_int field_size = 10;
//  gr_int grid_rank = 3;
//  // If grid rank is less than 3, set the other dimensions, 
//  // start indices, and end indices to 0.
//  gr_int grid_dimension[3], grid_start[3], grid_end[3];
//  for (int i = 0;i < 3;i++) {
//    grid_dimension[i] = 0; // the active dimension not including ghost zones.
//    grid_start[i] = 0;
//    grid_end[i] = 0;
//  }
//  grid_dimension[0] = field_size;
//  grid_end[0] = field_size - 1;
//
//  density       = new gr_float[field_size];
//  energy        = new gr_float[field_size];
//  x_velocity    = new gr_float[field_size];
//  y_velocity    = new gr_float[field_size];
//  z_velocity    = new gr_float[field_size];
//  // for primordial_chemistry >= 1
//  HI_density    = new gr_float[field_size];
//  HII_density   = new gr_float[field_size];
//  HeI_density   = new gr_float[field_size];
//  HeII_density  = new gr_float[field_size];
//  HeIII_density = new gr_float[field_size];
//  e_density     = new gr_float[field_size];
//  // for primordial_chemistry >= 2
//  HM_density    = new gr_float[field_size];
//  H2I_density   = new gr_float[field_size];
//  H2II_density  = new gr_float[field_size];
//  // for primordial_chemistry >= 3
//  DI_density    = new gr_float[field_size];
//  DII_density   = new gr_float[field_size];
//  HDI_density   = new gr_float[field_size];
//  // for metal_cooling = 1
//  metal_density = new gr_float[field_size];
//
//  // set temperature units
//  gr_float temperature_units =  mh * pow(units.a_units * 
//                                         units.length_units /
//                                         units.time_units, 2) / kboltz;
//
//  int i;
//  for (i = 0;i < field_size;i++) {
//    density[i] = 1.0;
//    HI_density[i] = chemistry_.HydrogenFractionByMass * density[i];
//    HII_density[i] = tiny_number * density[i];
//    HM_density[i] = tiny_number * density[i];
//    HeI_density[i] = (1.0 - chemistry_.HydrogenFractionByMass) * density[i];
//    HeII_density[i] = tiny_number * density[i];
//    HeIII_density[i] = tiny_number * density[i];
//    H2I_density[i] = tiny_number * density[i];
//    H2II_density[i] = tiny_number * density[i];
//    DI_density[i] = 2.0 * 3.4e-5 * density[i];
//    DII_density[i] = tiny_number * density[i];
//    HDI_density[i] = tiny_number * density[i];
//    e_density[i] = tiny_number * density[i];
//    // solar metallicity
//    metal_density[i] = chemistry_.SolarMetalFractionByMass * density[i];

//    x_velocity[i] = 0.0;
//    y_velocity[i] = 0.0;
//    z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
//    energy[i] = 1000. / temperature_units;
//  }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  // Evolving the chemistry_.
  // some timestep
//  gr_float dt = 3.15e7 * 1e6 / units.time_units;
//
//  if (solve_chemistry(chemistry, units,
//                      a_value, dt,
//                      grid_rank, grid_dimension,
//                      grid_start, grid_end,
//                      density, energy,
//                      x_velocity, y_velocity, z_velocity,
//                      HI_density, HII_density, HM_density,
//                      HeI_density, HeII_density, HeIII_density,
//                      H2I_density, H2II_density,
//                      DI_density, DII_density, HDI_density,
//                      e_density, metal_density) == 0) {
//    fprintf(stderr, "Error in solve_chemistry_.\n");
//    return 0;
//  }

// #include <string.h>

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodGrackle::EnzoMethodGrackle (EnzoConfig * c)
  : Method(), chemistry_()
{
  printf ("EnzoMethodGrackle()\n");

  /// Initialize parameters

  chemistry_.use_grackle            = true;
  chemistry_.Gamma                  = c->method_grackle_gamma;

  chemistry_.with_radiative_cooling = c->method_grackle_with_radiative_cooling;
  chemistry_.primordial_chemistry   = c->method_grackle_primordial_chemistry;
  chemistry_.metal_cooling          = c->method_grackle_metal_cooling;
  chemistry_.h2_on_dust             = c->method_grackle_h2_formation_on_dust;
  chemistry_.cmb_temperature_floor  = c->method_grackle_cmb_temperature_floor;
  chemistry_.grackle_data_file      = strdup(c->method_grackle_data_file.c_str());
  chemistry_.three_body_rate        = c->method_grackle_three_body_rate;
  chemistry_.cie_cooling            = c->method_grackle_cie_cooling;
  chemistry_.h2_optical_depth_approximation = c->method_grackle_h2_optical_depth_approximation;

  chemistry_.photoelectric_heating  = c->method_grackle_photoelectric_heating;

  chemistry_.photoelectric_heating_rate     
    = c->method_grackle_photoelectric_heating_rate;

  chemistry_.UVbackground           = c->method_grackle_UVbackground;

  chemistry_.UVbackground_table.Nz     = 0;
  chemistry_.UVbackground_table.z      = NULL;
  chemistry_.UVbackground_table.k24    = NULL;
  chemistry_.UVbackground_table.k25    = NULL;
  chemistry_.UVbackground_table.k26    = NULL;
  chemistry_.UVbackground_table.k27    = NULL;
  chemistry_.UVbackground_table.k28    = NULL;
  chemistry_.UVbackground_table.k29    = NULL;
  chemistry_.UVbackground_table.k30    = NULL;
  chemistry_.UVbackground_table.k31    = NULL;
  chemistry_.UVbackground_table.piHI   = NULL;
  chemistry_.UVbackground_table.piHeII = NULL;
  chemistry_.UVbackground_table.piHeI  = NULL;

  chemistry_.UVbackground_redshift_on   = c->method_grackle_UVbackground_redshift_on;
  chemistry_.UVbackground_redshift_off  = c->method_grackle_UVbackground_redshift_off;
  chemistry_.UVbackground_redshift_fullon  
    = c->method_grackle_UVbackground_redshift_fullon;
  chemistry_.UVbackground_redshift_drop = c->method_grackle_UVbackground_redshift_drop;

  chemistry_.Compton_xray_heating   = c->method_grackle_Compton_xray_heating;

  chemistry_.LWbackground_intensity = c->method_grackle_LWbackground_intensity;

  chemistry_.LWbackground_sawtooth_suppression 
    = c->method_grackle_LWbackground_sawtooth_suppression;

  chemistry_.HydrogenFractionByMass = c->method_grackle_HydrogenFractionByMass;
  chemistry_.DeuteriumToHydrogenRatio
    = c->method_grackle_DeuteriumToHydrogenRatio;
  chemistry_.SolarMetalFractionByMass     
    = c->method_grackle_SolarMetalFractionByMass;
  chemistry_.NumberOfTemperatureBins      
    = c->method_grackle_NumberOfTemperatureBins;
  chemistry_.ih2co                  = c->method_grackle_ih2co;
  chemistry_.ipiht                  = c->method_grackle_ipiht;
  chemistry_.TemperatureStart       = c->method_grackle_TemperatureStart;
  chemistry_.TemperatureEnd         = c->method_grackle_TemperatureEnd;
  chemistry_.comp_xray              = c->method_grackle_comp_xray;
  chemistry_.temp_xray              = c->method_grackle_temp_xray;
  chemistry_.CaseBRecombination     = c->method_grackle_CaseBRecombination;
  chemistry_.NumberOfDustTemperatureBins  
    = c->method_grackle_NumberOfDustTemperatureBins;
  chemistry_.DustTemperatureStart   = c->method_grackle_DustTemperatureStart;
  chemistry_.DustTemperatureEnd     = c->method_grackle_DustTemperatureEnd;

  chemistry_.cloudy_metal.grid_rank = c->mesh_root_rank;
  chemistry_.cloudy_electron_fraction_factor 
    = c->method_grackle_cloudy_electron_fraction_factor;
}

//----------------------------------------------------------------------

void EnzoMethodGrackle::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute ( CommBlock * comm_block) throw()
{
  initialize_(comm_block);

  /// Initialize Fields
//  gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity,
//    *HI_density, *HII_density, *HM_density,
//    *HeI_density, *HeII_density, *HeIII_density,
//    *H2I_density, *H2II_density,
//    *DI_density, *DII_density, *HDI_density,
//    *e_density, *metal_density;



}

//----------------------------------------------------------------------

double EnzoMethodGrackle::timestep ( CommBlock * comm_block ) throw()
{
  initialize_(comm_block);
  return std::numeric_limits<double>::max();
}

//======================================================================
