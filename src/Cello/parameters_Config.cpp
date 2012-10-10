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
  p | field_refresh_faces;

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

  num_fields = parameters->list_length("Field:fields"); 

  field_fields.resize(num_fields);
  for (int i=0; i<num_fields; i++) {
    field_fields[i] = parameters->list_value_string(i, "Field:fields");
  }

  if (parameters->type("Field:ghosts") == parameter_integer) {
    field_ghosts[0] = parameters->value_integer("Field:ghosts",0);
    field_ghosts[1] = parameters->value_integer("Field:ghosts",0);
    field_ghosts[2] = parameters->value_integer("Field:ghosts",0);
  } else if (parameters->type("Field:ghosts") == parameter_list) {
    field_ghosts[0] = parameters->list_value_integer(0,"Field:ghosts",0);
    field_ghosts[1] = parameters->list_value_integer(1,"Field:ghosts",0);
    field_ghosts[2] = parameters->list_value_integer(2,"Field:ghosts",0);
  }

  field_alignment = parameters->value_integer("Field:alignment",8);

  field_centering[0].resize(num_fields);
  field_centering[1].resize(num_fields);
  field_centering[2].resize(num_fields);
  for (int i=0; i<num_fields; i++) {

    std::string param_name = 
      std::string("Field:") + field_fields[i] + ":centering";

    field_centering[0][i] = parameters->list_value_logical(0,param_name,true);
    field_centering[1][i] = parameters->list_value_logical(1,param_name,true);
    field_centering[2][i] = parameters->list_value_logical(2,param_name,true);
    
  }

  field_courant = parameters->value_float  ("Field:courant",0.6);

  field_padding = parameters->value_integer("Field:padding",0);

  // Field precision

  std::string precision_str = parameters->value_string("Field:precision","default");

  if      (precision_str == "default")   field_precision = precision_default;
  else if (precision_str == "single")    field_precision = precision_single;
  else if (precision_str == "double")    field_precision = precision_double;
  else if (precision_str == "quadruple") field_precision = precision_quadruple;
  else {
    ERROR1 ("Simulation::initialize_data_descr_()", 
	    "Unknown precision %s",
	    precision_str.c_str());
  }

  field_refresh_corners = parameters->value_logical ("Field:refresh:corners",true);
  field_refresh_edges   = parameters->value_logical ("Field:refresh:edges",  true);
  field_refresh_faces   = parameters->value_logical ("Field:refresh:faces",  true);

  mesh_root_rank = parameters->value_integer("Mesh:root_rank",0);

  if (mesh_root_rank < 2) field_ghosts[1] = 0;
  if (mesh_root_rank < 3) field_ghosts[2] = 0;
  

  initial_cycle = parameters->value_integer("Initial:cycle",0);
  initial_time  = parameters->value_float  ("Initial:time",0.0);
  initial_type  = parameters->value_string("Initial:type","default");

  //  initial_name;
  //  initial_value

  mesh_root_blocks[0] = parameters->list_value_integer(0,"Mesh:root_blocks",1);
  mesh_root_blocks[1] = parameters->list_value_integer(1,"Mesh:root_blocks",1);
  mesh_root_blocks[2] = parameters->list_value_integer(2,"Mesh:root_blocks",1);

#ifndef CONFIG_USE_CHARM
  int root_blocks = mesh_root_blocks[0]*mesh_root_blocks[1]*mesh_root_blocks[2];
  GroupProcess * group_process = GroupProcess::create();
  ASSERT4 ("Simulation::initialize_hierarchy_",
	   "Product of Mesh:root_blocks = [%d %d %d] must equal MPI_Comm_size",
	   mesh_root_blocks[0],
	   mesh_root_blocks[1],
	   mesh_root_blocks[2],
	   group_process->size(),
	   root_blocks==group_process->size());
  delete group_process;
#endif


  mesh_root_size[0] = parameters->list_value_integer(0,"Mesh:root_size",1);
  mesh_root_size[1] = parameters->list_value_integer(1,"Mesh:root_size",1);
  mesh_root_size[2] = parameters->list_value_integer(2,"Mesh:root_size",1);

  //  method_sequence;

  monitor_debug = parameters->value_logical("Monitor:debug",false);

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

  stopping_cycle = parameters->value_integer
    ( "Stopping:cycle" , std::numeric_limits<int>::max() );
  stopping_time  = parameters->value_float
    ( "Stopping:time" , std::numeric_limits<double>::max() );

  //  testing_cycle_final;
  //  testing_time_final;

  timestep_type = parameters->value_string("Timestep:type","default");

}

//======================================================================

