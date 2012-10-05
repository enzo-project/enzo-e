// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-03
/// @brief    Implementation of the Config class 

#include "cello.hpp"
#include "parameters.hpp"

//----------------------------------------------------------------------

Config::Config() throw ()
{
}

//----------------------------------------------------------------------

Config::~Config() throw ()
{
}

//----------------------------------------------------------------------

Config::Config(const Config & config) throw ()
/// @param     config  Object being copied
{
}

//----------------------------------------------------------------------

Config & Config::operator= (const Config & config) throw ()
/// @param     config  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE("Config::operator=");
  return *this;
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Config::pup (PUP::er &p)
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

  p | field_alignment;
  PUParray(p,field_centering,MAX_FIELDS);
  p | field_courant;
  p | field_fields;
  PUParray(p,field_ghosts,3);
  p | field_padding;
  p | field_precision;
  p | field_refresh_corners;
  p | field_refresh_edges;

  p | initial_cycle;
  p | initial_type;
  p | initial_time;
  p | initial_name;
  PUParray(p,initial_value,MAX_FIELDS);

  PUParray(p,mesh_root_blocks,3);
  p | mesh_root_rank;
  PUParray(p,mesh_root_size,3);

  p | method_sequence;

  p | monitor_debug;

  p | output_file_groups;
  PUParray (p,output_axis,MAX_FILE_GROUPS);
  PUParray (p,output_colormap,MAX_FILE_GROUPS);
  PUParray (p,output_colormap_alpha,MAX_FILE_GROUPS);
  PUParray (p,output_field_list,MAX_FILE_GROUPS);
  PUParray (p,output_name,MAX_FILE_GROUPS);
  PUParray (p,output_dir,MAX_FILE_GROUPS);
  PUParray (p,output_schedule,MAX_FILE_GROUPS);
  PUParray (p,output_type,MAX_FILE_GROUPS);

  p | physics_cosmology;
  p | physics_cosmology_comoving_box_size;
  p | physics_cosmology_hubble_constant_now;
  p | physics_cosmology_initial_redshift;
  p | physics_cosmology_max_expansion_rate;
  p | physics_cosmology_omega_lamda_now;
  p | physics_cosmology_omega_matter_now;
  p | physics_gamma;

  p | stopping_cycle;
  p | stopping_time;

  p | testing_cycle_final;
  p | testing_time_final;

  p | timestep_type;

}

#endif

//----------------------------------------------------------------------

void Config::read(Parameters * parameters) throw()
{
  boundary_type = parameters->value_string("Boundary:type","");

  for (int i=0; i<3; i++)  {
    domain_lower[i] = parameters->list_value_float(i, "Domain:lower", 0.0);
    domain_upper[i] = parameters->list_value_float(i, "Domain:upper", 0.0);
  }

  //  enzo_ppm_density_floor;
  //  enzo_ppm_diffusion;
  //  enzo_ppm_dual_energy;
  //  enzo_ppm_dual_energy_eta_1;
  //  enzo_ppm_dual_energy_eta_2;
  //  enzo_ppm_flattening;
  //  enzo_ppm_minimum_pressure_support_parameter;
  //  enzo_ppm_number_density_floor;
  //  enzo_ppm_pressure_floor;
  //  enzo_ppm_pressure_free;
  //  enzo_ppm_steepening;
  //  enzo_ppm_temperature_floor;
  //  enzo_ppm_use_minimum_pressure_support;

  //  field_alignment;
  //  field_centering
  //  field_courant;
  //  field_fields;
  //  field_ghosts
  //  field_padding;
  //  field_precision;
  //  field_refresh_corners;
  //  field_refresh_edges;

  initial_cycle = parameters->value_integer("Initial:cycle",0);
  initial_time  = parameters->value_float  ("Initial:time",0.0);
  //  initial_type;
  //  initial_name;
  //  initial_value

  //  mesh_root_blocks
  mesh_root_rank = parameters->value_integer("Mesh:root_rank",0);
  //  mesh_root_size

  //  method_sequence;

  //  monitor_debug;

  //  output_file_groups;
  //  output_axis
  //  output_colormap
  //  output_colormap_alpha
  //  output_field_list
  //  output_name
  //  output_dir
  //  output_schedule
  //  output_type

  // RENAME physics_ as enzo_

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

  //  stopping_cycle;
  //  stopping_time;

  //  testing_cycle_final;
  //  testing_time_final;

  //  timestep_type;

}

//======================================================================

