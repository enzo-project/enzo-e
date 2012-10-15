// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_EnzoConfig.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-03
/// @brief    Implementation of the EnzoConfig class 

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoConfig::EnzoConfig() throw ()
{
}

//----------------------------------------------------------------------

EnzoConfig::~EnzoConfig() throw ()
{
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void EnzoConfig::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
  p | boundary_type;

  PUParray(p,domain_lower,3);
  PUParray(p,domain_upper,3);

  p | enzo_ppm_density_floor;
  p | enzo_ppm_diffusion;
  p | enzo_ppm_dual_energy;
  p | enzo_ppm_dual_energy_eta_1;
  p | enzo_ppm_dual_energy_eta_2;
  p | enzo_ppm_flattening;
  p | enzo_ppm_minimum_pressure_support_parameter;
  p | enzo_ppm_number_density_floor;
  p | enzo_ppm_pressure_floor;
  p | enzo_ppm_pressure_free;
  p | enzo_ppm_steepening;
  p | enzo_ppm_temperature_floor;
  p | enzo_ppm_use_minimum_pressure_support;

  p | physics_cosmology;
  p | physics_cosmology_comoving_box_size;
  p | physics_cosmology_hubble_constant_now;
  p | physics_cosmology_initial_redshift;
  p | physics_cosmology_max_expansion_rate;
  p | physics_cosmology_omega_lamda_now;
  p | physics_cosmology_omega_matter_now;
  p | physics_gamma;

}

#endif

//----------------------------------------------------------------------

void EnzoConfig::read(Parameters * parameters) throw()
{

  TRACE("BEGIN EnzoConfig::read()");

  // Read Cello parameters

  ((Config*)this) -> read (parameters);

  double floor_default = 1e-6;

  enzo_ppm_density_floor = parameters->value_float
    ("Enzo:ppm:density_floor",  floor_default);
  enzo_ppm_diffusion = parameters->value_logical 
    ("Enzo:ppm:diffusion",  false);
  enzo_ppm_dual_energy = parameters->value_logical 
    ("Enzo:ppm:dual_energy",false);
  enzo_ppm_dual_energy_eta_1 = parameters->value_float
    ("Enzo:ppm:dual_energy_eta_1", 0.001);
  enzo_ppm_dual_energy_eta_2 = parameters->value_float
    ("Enzo:ppm:dual_energy_eta_2", 0.1);
  enzo_ppm_flattening = parameters->value_integer
    ("Enzo:ppm:flattening", 3);
  enzo_ppm_minimum_pressure_support_parameter = parameters->value_integer
    ("Enzo:ppm:minimum_pressure_support_parameter",100);
  enzo_ppm_number_density_floor = parameters->value_float
    ("Enzo:ppm:number_density_floor", floor_default);
  enzo_ppm_pressure_floor = parameters->value_float
    ("Enzo:ppm:pressure_floor", floor_default);
  enzo_ppm_pressure_free = parameters->value_logical
    ("Enzo:ppm:pressure_free",false);
  enzo_ppm_steepening = parameters->value_logical 
    ("Enzo:ppm:steepening", false);
  enzo_ppm_temperature_floor = parameters->value_float
    ("Enzo:ppm:temperature_floor", floor_default);
  enzo_ppm_use_minimum_pressure_support = parameters->value_logical
    ("Enzo:ppm:use_minimum_pressure_support",false);


  physics_cosmology = parameters->value_logical ("Physics:cosmology",false);
  physics_cosmology_comoving_box_size = parameters->value_float
    ("Physics:cosmology:comoving_box_size", 64.0);
  physics_cosmology_hubble_constant_now = parameters->value_float
    ("Physics:cosmology:hubble_constant_now",0.701);
  physics_cosmology_initial_redshift = parameters->value_float
    ("Physics:cosmology:initial_redshift",  20.0);;
  physics_cosmology_max_expansion_rate = parameters->value_float
    ("Physics:cosmology:max_expansion_rate", 0.01);
  physics_cosmology_omega_lamda_now = parameters->value_float
    ("Physics:cosmology:omega_lambda_now",   0.721);
  physics_cosmology_omega_matter_now = parameters->value_float
    ("Physics:cosmology:omega_matter_now",   0.279);
  physics_gamma = parameters->value_float ("Physics:gamma",5.0/3.0);

  TRACE("END   EnzoConfig::read()");
}

//======================================================================

