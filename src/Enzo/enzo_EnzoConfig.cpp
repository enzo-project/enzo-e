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

void EnzoConfig::pup (PUP::er &p)
{
  Config::pup(p);
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  p | ppm_density_floor;
  p | ppm_diffusion;
  p | ppm_dual_energy;
  p | ppm_dual_energy_eta_1;
  p | ppm_dual_energy_eta_2;
  p | ppm_flattening;
  p | ppm_minimum_pressure_support_parameter;
  p | ppm_number_density_floor;
  p | ppm_pressure_floor;
  p | ppm_pressure_free;
  p | ppm_steepening;
  p | ppm_temperature_floor;
  p | ppm_use_minimum_pressure_support;

  p | field_gamma;

  p | cosmology;
  p | cosmology_comoving_box_size;
  p | cosmology_hubble_constant_now;
  p | cosmology_initial_redshift;
  p | cosmology_max_expansion_rate;
  p | cosmology_omega_lamda_now;
  p | cosmology_omega_matter_now;

  PUParray(p,sedov_array,3);
  p | sedov_radius_relative;
  p | sedov_pressure_in;
  p | sedov_pressure_out;
  p | sedov_density;

  p | interpolation_method;

  p | method_heat_alpha;

}

//----------------------------------------------------------------------

void EnzoConfig::read(Parameters * parameters) throw()
{

  TRACE("BEGIN EnzoConfig::read()");

  // Read Cello parameters

  
  TRACE("EnzoCharm::read calling Config::read()");
  ((Config*)this) -> read (parameters);

  double floor_default = 1e-6;

  ppm_density_floor = parameters->value_float
    ("Method:ppm:density_floor",  floor_default);
  ppm_diffusion = parameters->value_logical 
    ("Method:ppm:diffusion",  false);
  ppm_dual_energy = parameters->value_logical 
    ("Method:ppm:dual_energy",false);
  ppm_dual_energy_eta_1 = parameters->value_float
    ("Method:ppm:dual_energy_eta_1", 0.001);
  ppm_dual_energy_eta_2 = parameters->value_float
    ("Method:ppm:dual_energy_eta_2", 0.1);
  ppm_flattening = parameters->value_integer
    ("Method:ppm:flattening", 3);
  ppm_minimum_pressure_support_parameter = parameters->value_integer
    ("Method:ppm:minimum_pressure_support_parameter",100);
  ppm_number_density_floor = parameters->value_float
    ("Method:ppm:number_density_floor", floor_default);
  ppm_pressure_floor = parameters->value_float
    ("Method:ppm:pressure_floor", floor_default);
  ppm_pressure_free = parameters->value_logical
    ("Method:ppm:pressure_free",false);
  ppm_steepening = parameters->value_logical 
    ("Method:ppm:steepening", false);
  ppm_temperature_floor = parameters->value_float
    ("Method:ppm:temperature_floor", floor_default);
  ppm_use_minimum_pressure_support = parameters->value_logical
    ("Method:ppm:use_minimum_pressure_support",false);


  cosmology = parameters->value_logical ("Method:cosmology",false);
  cosmology_comoving_box_size = parameters->value_float
    ("Method:cosmology:comoving_box_size", 64.0);
  cosmology_hubble_constant_now = parameters->value_float
    ("Method:cosmology:hubble_constant_now",0.701);
  cosmology_initial_redshift = parameters->value_float
    ("Method:cosmology:initial_redshift",  20.0);;
  cosmology_max_expansion_rate = parameters->value_float
    ("Method:cosmology:max_expansion_rate", 0.01);
  cosmology_omega_lamda_now = parameters->value_float
    ("Method:cosmology:omega_lambda_now",   0.721);
  cosmology_omega_matter_now = parameters->value_float
    ("Method:cosmology:omega_matter_now",   0.279);

  field_gamma = parameters->value_float ("Field:gamma",5.0/3.0);

  TRACE1("field_gamma = %f",field_gamma);

  sedov_array[0] = parameters->list_value_integer (0,"Initial:sedov:array",1);
  sedov_array[1] = parameters->list_value_integer (1,"Initial:sedov:array",1);
  sedov_array[2] = parameters->list_value_integer (2,"Initial:sedov:array",1);

  sedov_radius_relative = 
    parameters->value_float("Initial:sedov:radius_relative",0.1);
  sedov_pressure_in = 
    parameters->value_float("Initial:sedov:pressure_in",1.0);
  sedov_pressure_out = 
    parameters->value_float("Initial:sedov:pressure_in",1e-5);
  sedov_density = 
    parameters->value_float("Initial:sedov:density",1.0);

  interpolation_method = parameters->value_string 
    ("Field:interpolation_method","SecondOrderA");


  method_heat_alpha = parameters->value_float 
    ("Method:heat:alpha",1.0);

  TRACE("END   EnzoConfig::read()");
}

//======================================================================

