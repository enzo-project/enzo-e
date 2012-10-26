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
  Config::pup(p);
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

  p | enzo_cosmology;
  p | enzo_cosmology_comoving_box_size;
  p | enzo_cosmology_hubble_constant_now;
  p | enzo_cosmology_initial_redshift;
  p | enzo_cosmology_max_expansion_rate;
  p | enzo_cosmology_omega_lamda_now;
  p | enzo_cosmology_omega_matter_now;
  p | enzo_gamma;

  PUParray(p,enzo_sedov_array,3);

}

#endif

//----------------------------------------------------------------------

void EnzoConfig::read(Parameters * parameters) throw()
{

  TRACE("BEGIN EnzoConfig::read()");

  // Read Cello parameters

  
  TRACE("EnzoCharm::read calling Config::read()");
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


  enzo_cosmology = parameters->value_logical ("Enzo:cosmology",false);
  enzo_cosmology_comoving_box_size = parameters->value_float
    ("Enzo:cosmology:comoving_box_size", 64.0);
  enzo_cosmology_hubble_constant_now = parameters->value_float
    ("Enzo:cosmology:hubble_constant_now",0.701);
  enzo_cosmology_initial_redshift = parameters->value_float
    ("Enzo:cosmology:initial_redshift",  20.0);;
  enzo_cosmology_max_expansion_rate = parameters->value_float
    ("Enzo:cosmology:max_expansion_rate", 0.01);
  enzo_cosmology_omega_lamda_now = parameters->value_float
    ("Enzo:cosmology:omega_lambda_now",   0.721);
  enzo_cosmology_omega_matter_now = parameters->value_float
    ("Enzo:cosmology:omega_matter_now",   0.279);
  enzo_gamma = parameters->value_float ("Enzo:gamma",5.0/3.0);

  TRACE1("enzo_gamma = %f",enzo_gamma);

  enzo_sedov_array[0] = parameters->list_value_integer (0,"Enzo:sedov:array",1);
  enzo_sedov_array[1] = parameters->list_value_integer (1,"Enzo:sedov:array",1);
  enzo_sedov_array[2] = parameters->list_value_integer (2,"Enzo:sedov:array",1);

  TRACE("END   EnzoConfig::read()");
}

//======================================================================

