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

  p | num_fields;
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

  p | num_file_groups;
  p | output_file_groups;
  PUParray (p,output_type,MAX_FILE_GROUPS);
  PUParray (p,output_image_axis,MAX_FILE_GROUPS);
  PUParray (p,output_image_colormap,MAX_FILE_GROUPS);
  PUParray (p,output_image_colormap_alpha,MAX_FILE_GROUPS);
  PUParray (p,output_field_list,MAX_FILE_GROUPS);
  PUParray (p,output_stride,MAX_FILE_GROUPS);
  PUParray (p,output_name,MAX_FILE_GROUPS);
  PUParray (p,output_dir,MAX_FILE_GROUPS);
  PUParray (p,output_schedule_type,MAX_FILE_GROUPS);
  PUParray (p,output_schedule_var,MAX_FILE_GROUPS);
  PUParray (p,output_schedule_start,MAX_FILE_GROUPS);
  PUParray (p,output_schedule_stop,MAX_FILE_GROUPS);
  PUParray (p,output_schedule_step,MAX_FILE_GROUPS);
  PUParray (p,output_schedule_list,MAX_FILE_GROUPS);

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

  ASSERT2 ("Config::read","Number of fields %d exceeds MAX_FIELDS %d",
	   num_fields, MAX_FIELDS, num_fields <= MAX_FIELDS);

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
    ERROR1 ("Config::read()", "Unknown precision %s",
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
  ASSERT4 ("Config::read()",
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

  int size = parameters->list_length("Method:sequence");
  method_sequence.resize(size);
  for (int i=0; i<size; i++) {
    method_sequence[i] = parameters->list_value_string(i,"Method:sequence");
  }

  monitor_debug = parameters->value_logical("Monitor:debug",false);

  // Output

  parameters->group_set(0,"Output");

  num_file_groups = parameters->list_length("file_groups");

  ASSERT2 ("Config::read","Number of file groups %d exceeds MAX_FILE_GROUPS %d",
	   num_file_groups, MAX_FILE_GROUPS, num_file_groups <= MAX_FILE_GROUPS);

  parameters->group_set(0,"Output");


  output_file_groups.resize(num_file_groups);

  for (int index=0; index<num_file_groups; index++) {

    output_file_groups[index] = 
      parameters->list_value_string (index,"Output:file_groups","unknown");

    parameters->group_set(1,output_file_groups[index]);

    output_type[index] = parameters->value_string("type","unknown");

    if (output_type[index] == "unknown") {
      ERROR1("Config::read",
	     "Output:%s:type parameter is undefined",
	     output_file_groups[index].c_str());
    }

    output_stride[index] = parameters->value_integer("stride",0);

    if (parameters->type("dir") == parameter_string) {
      output_dir[index].resize(1);
      output_dir[index][0] = parameters->value_string("dir","");
    } else if (parameters->type("dir") == parameter_list) {
      int size = parameters->list_length("dir");
      if (size > 0) output_dir[index].resize(size);
      for (int i=0; i<size; i++) {
	output_dir[index][i] = parameters->list_value_string(i,"dir","");
      }
    }

    if (parameters->type("name") == parameter_string) {
      output_name[index].resize(1);
      output_name[index][0] = parameters->value_string("name","");
    } else if (parameters->type("name") == parameter_list) {
      int size = parameters->list_length("name");
      if (size > 0) output_name[index].resize(size);
      for (int i=0; i<size; i++) {
	output_name[index][i] = parameters->list_value_string(i,"name","");
      }
    }
    //  output_axis
    //  output_colormap
    //  output_colormap_alpha
    if (parameters->type("field_list") == parameter_list) {
      int length = parameters->list_length("field_list");
      output_field_list[index].resize(length);
      for (int i=0; i<length; i++) {
	output_field_list[index][i] = parameters->list_value_string(i,"field_list","");
      }
    }
      
    //  output_schedule

    ASSERT1("Config::read",
	   "'schedule' is not defined for output file group %s",
	    output_file_groups[index].c_str(),
	    (parameters->type("schedule") != parameter_unknown));

    int length = parameters->list_length("schedule");
    ASSERT1("Config::read","Incorrect 'schedule' for Output group %s",
	    output_file_groups[index].c_str(), (length >= 3));

    output_schedule_var[index]  = parameters->list_value_string(0,"schedule");
    output_schedule_type[index] = parameters->list_value_string(1,"schedule");

    bool var_is_int;

    if (output_schedule_var[index] == "cycle") {
      var_is_int = true;
    } else if (output_schedule_var[index] == "time") {
      var_is_int = false;
    } else {
      ERROR2 ("Config::read",
	      "Schedule variable %s is not recognized for output file group %s",
	      output_schedule_var[index].c_str(),output_file_groups[index].c_str());
    }

    const int    max_int    = std::numeric_limits<int>::max();
    const double max_double = std::numeric_limits<double>::max();

    if (output_schedule_type[index] == "interval") {
      if (length == 3) { // step
	if (var_is_int) {
	  output_schedule_start[index] = 0;
	  output_schedule_step[index]  = parameters->list_value_integer(2,"schedule");
	  output_schedule_stop[index]  = max_int;
	} else {
	  output_schedule_start[index] = 0.0;
	  output_schedule_step[index]  = parameters->list_value_float(2,"schedule");
	  output_schedule_stop[index]  = max_double;
	}
      } else if (length == 4) { // start, step
	if (var_is_int) {
	  output_schedule_start[index] = parameters->list_value_integer(2,"schedule");
	  output_schedule_step[index]  = parameters->list_value_integer(3,"schedule");
	  output_schedule_stop[index]  = max_int;
	} else {
	  output_schedule_start[index] = parameters->list_value_float(2,"schedule");
	  output_schedule_step[index]  = parameters->list_value_float(3,"schedule");
	  output_schedule_stop[index]  = max_double;
	}
      } else if (length == 5) { // start, step, stop
	if (var_is_int) {
	  output_schedule_start[index] = parameters->list_value_integer(2,"schedule");
	  output_schedule_step[index]  = parameters->list_value_integer(3,"schedule");
	  output_schedule_stop[index]  = parameters->list_value_integer(4,"schedule");
	} else {
	  output_schedule_start[index] = parameters->list_value_float(2,"schedule");
	  output_schedule_step[index]  = parameters->list_value_float(3,"schedule");
	  output_schedule_stop[index]  = parameters->list_value_float(4,"schedule");
	}
      }
    } else if (output_schedule_type[index] == "list") {
      output_schedule_list[index].resize(length-2);
      for (int i=2; i<length; i++) {
	if (var_is_int) {
	  output_schedule_list[index][i-2] = parameters->list_value_integer(i,"schedule");
	} else {
	  output_schedule_list[index][i-2] = parameters->list_value_float(i,"schedule");
	}
      }
    } else {
      ERROR2 ("Config::read",
	      "Schedule type %s is not recognized for output file group %s",
	      output_schedule_type[index].c_str(),output_file_groups[index].c_str());
    }

    // Image 

    TRACE2 ("output_type[%d] = %s",index,output_type[index].c_str());

    if (output_type[index] == "image") {

      output_image_axis[index] = "z";

      if (parameters->type("axis") != parameter_unknown) {
	std::string axis = parameters->value_string("axis");
	ASSERT2("Problem::initialize_output",
		"Output %s axis %d must be \"x\", \"y\", or \"z\"",
		output_file_groups[index].c_str(), axis.c_str(),
		axis=="x" || axis=="y" || axis=="z");
      } 

      if (parameters->type("colormap") == parameter_list) {
	int size = parameters->list_length("colormap");
	output_image_colormap[index].resize(size);
	for (int i=0; i<size; i++) {
	  double cm = 
	    output_image_colormap[index][i] = parameters->list_value_float(i,"colormap",0.0);
	  printf ("cm[%d] = %f\n",i,cm);
	}
      }

      if (parameters->type("colormap_alpha") == parameter_list) {
	int size = parameters->list_length("colormap_alpha");
	output_image_colormap_alpha[index].resize(size);
	for (int i=0; i<size; i++) {
	  output_image_colormap_alpha[index][i] = 
	    parameters->list_value_float(i,"colormap_alpha",0.0);
	}
      }

    }

  }  
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

