// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Config.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-03
/// @brief    Implementation of the Config class 
///
/// Last review of parameters was 2015-09-10 with revision -r 3836 

#include "cello.hpp"
#include "parameters.hpp"
#include <iostream>

Config g_config;

//----------------------------------------------------------------------

void Config::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  PUP::able::pup(p);

  TRACEPUP;

  // Adapt

  p | num_adapt;
  p | adapt_list;
  p | adapt_interval;
  p | adapt_min_face_rank;
  p | adapt_type;
  p | adapt_field_list;
  p | adapt_min_refine;
  p | adapt_max_coarsen;
  p | adapt_min_refine2;
  p | adapt_max_coarsen2;
  p | adapt_max_level;
  p | adapt_level_exponent;
  p | adapt_include_ghosts;
  p | adapt_output;
  p | adapt_schedule_index;
  
  // Balance

  p | balance_schedule_index;
  p | balance_type;

  // Boundary

  p | num_boundary;
  p | boundary_list;
  p | boundary_type;
  p | boundary_axis;
  p | boundary_face;
  p | boundary_mask;
  p | boundary_field_list;

  // Domain

  PUParray(p,domain_lower,3);
  PUParray(p,domain_upper,3);

  // Field

  p | num_fields;
  p | field_list;
  p | field_alignment;
  PUParray(p,field_centering,3);
  PUParray(p,field_ghost_depth,3);
  p | field_padding;
  p | field_history;
  p | field_precision;
  p | field_prolong;
  p | field_restrict;
  p | field_group_list;

  // Initial

  p | num_initial;
  p | initial_new;
  p | initial_list;
  p | initial_cycle;
  p | initial_time;

  p | initial_restart;
  p | initial_restart_dir;

  p | initial_trace_name;
  p | initial_trace_field;
  p | initial_trace_mpp;
  p | initial_trace_dx;
  p | initial_trace_dy;
  p | initial_trace_dz;

  // Memory

  p | memory_active;
  p | memory_warning_mb;
  p | memory_limit_gb;

  // Mesh

  PUParray(p,mesh_root_blocks,3);
  p | mesh_root_rank;
  PUParray(p,mesh_root_size,3);
  p | mesh_min_level;
  p | mesh_max_level;
  p | mesh_max_initial_level;

  // Method

  p | num_method;
  p | method_courant_global;
  p | method_list;
  p | method_schedule_index;
  p | method_file_name;
  p | method_path_name;
  p | method_close_files_seconds_stagger;
  p | method_close_files_seconds_delay;
  p | method_close_files_group_size;
  p | method_courant;
  p | method_debug_print;
  p | method_debug_coarse;
  p | method_debug_ghost;
  p | method_flux_correct_group;
  p | method_flux_correct_enable;
  p | method_flux_correct_min_digits_fields;
  p | method_flux_correct_min_digits_values;
  p | method_flux_correct_single_array;
  p | method_field_list;
  p | method_particle_list;
  PUParray (p,method_output_blocking,3);
  p | method_output_all_blocks;
  p | method_prolong;
  p | method_ghost_depth;
  p | method_min_face_rank;
  p | method_all_fields;
  p | method_all_particles;

  p | method_timestep;
  p | method_trace_name;
  p | method_type;
  p | method_null_dt;

  // Monitor

  p | monitor_debug;
  p | monitor_verbose;

  // Output

  p | num_output;
  p | output_list;
  p | output_type;
  p | output_axis;
  p | output_colormap;
  p | output_image_lower;
  p | output_image_upper;
  p | output_image_type;
  p | output_image_log;
  p | output_image_abs;
  p | output_image_mesh_color;
  p | output_image_mesh_order;
  p | output_image_color_particle_attribute;
  p | output_image_size;
  p | output_image_reduce_type;
  p | output_image_ghost;
  p | output_image_face_rank;
  p | output_image_min;
  p | output_image_max;
  p | output_min_level;
  p | output_max_level;
  p | output_leaf_only;
  p | output_schedule_index;
  p | output_dir;
  p | output_dir_global;
  p | output_stride_write;
  p | output_stride_wait;
  p | output_field_list;
  p | output_particle_list;
  p | output_checkpoint_file;
  p | output_name;
  p | index_schedule;
  p | schedule_list;
  p | schedule_type;
  p | schedule_var;
  p | schedule_start;
  p | schedule_stop;
  p | schedule_step;

  // Particles

  p | num_particles;
  p | particle_list;
  p | particle_index;
  p | particle_interleaved;
  p | particle_constant_name;
  p | particle_constant_type;
  p | particle_constant_value;
  p | particle_attribute_name;
  p | particle_attribute_type;
  PUParray (p,particle_attribute_position,3);
  PUParray (p,particle_attribute_velocity,3);
  p | particle_batch_size;
  p | particle_group_list;

  // Performance

  p | performance_papi_counters;
  p | performance_projections_on_at_start;
  p | performance_warnings;
  p | performance_on_schedule_index;
  p | performance_off_schedule_index;

  // Physics
  
  p | num_physics;
  p | physics_list;

  // Solvers
  
  p | num_solvers;
  p | solver_list;
  p | solver_index;
  p | solver_type;
  p | solver_solve_type;
  p | solver_iter_max;
  p | solver_res_tol;
  p | solver_diag_precon;
  p | solver_monitor_iter;
  p | solver_restrict;
  p | solver_prolong;
  p | solver_min_level;
  p | solver_max_level;
  p | solver_field_x;
  p | solver_field_b;
  
  // Stopping

  p | stopping_cycle;
  p | stopping_time;
  p | stopping_seconds;
  p | stopping_interval;

  // Testing

  p | testing_cycle_final;
  p | testing_time_final;
  p | testing_time_tolerance;

}

//----------------------------------------------------------------------

void Config::read(Parameters * p) throw()
{
  TRACE("BEGIN Config::read()");
  read_adapt_(p);
  read_balance_(p);
  read_boundary_(p);
  read_domain_(p);
  read_field_(p);
  read_initial_(p);
  read_memory_(p);
  read_mesh_(p);
  read_solver_(p); // before read_method_
  read_method_(p);
  read_monitor_(p);
  read_output_(p);
  read_particle_(p);
  read_performance_(p);
  read_physics_(p);
  read_stopping_(p);
  read_testing_(p);
  read_units_(p);

  TRACE("END   Config::read()");

}

//----------------------------------------------------------------------

void Config::read_adapt_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Adapt
  //--------------------------------------------------

  adapt_interval = p->value ("Adapt:interval",1);

  //--------------------------------------------------

  num_adapt = p->list_length("Adapt:list");

  adapt_list           .resize(num_adapt);
  adapt_type           .resize(num_adapt);
  adapt_field_list     .resize(num_adapt);
  adapt_min_refine     .resize(num_adapt);
  adapt_max_coarsen    .resize(num_adapt);
  adapt_min_refine2    .resize(num_adapt);
  adapt_max_coarsen2   .resize(num_adapt);
  adapt_max_level      .resize(num_adapt);
  adapt_level_exponent .resize(num_adapt);
  adapt_include_ghosts .resize(num_adapt);
  adapt_output         .resize(num_adapt);
  adapt_schedule_index .resize(num_adapt);

  adapt_min_face_rank = p->value_integer("Adapt:min_face_rank",0);

  for (int ia=0; ia<num_adapt; ia++) {

    adapt_list[ia] = p->list_value_string (ia,"Adapt:list","unknown");

    std::string prefix = "Adapt:" + adapt_list[ia] + ":";

    adapt_type[ia] = p->value_string(prefix+"type","unknown");

    std::string param_str = prefix + "field_list";

    int type = p->type(param_str);

    if (type == parameter_list) {
      const int n = p->list_length(param_str);
      adapt_field_list[ia].resize(n);
      for (int index=0; index<n; index++) {
	adapt_field_list[ia][index] = p->value(index,param_str,"none");
      }
    } else if (type == parameter_string) {
      adapt_field_list[ia].resize(1);
      adapt_field_list[ia][0] = p->value_string (param_str,"none");
    } else if (type != parameter_unknown) {
      ERROR2 ("Config::read()", "Incorrect parameter type %d for %s",
	      type,param_str.c_str());
    }

    //--------------------------------------------------

    if (p->type(prefix + "min_refine") == parameter_list) {

      adapt_min_refine[ia]  = p->list_value_float (0,prefix + "min_refine",0.3);
      adapt_max_coarsen[ia] = p->list_value_float (0,prefix + "max_coarsen",
						  0.5*adapt_min_refine[ia]);

      adapt_min_refine2[ia]  = p->list_value_float (0,prefix + "min_refine",0.3);
      adapt_max_coarsen2[ia] = p->list_value_float (1,prefix + "max_coarsen",
						   0.5*adapt_min_refine2[ia]);
    } else {

      adapt_min_refine[ia]  = p->value (prefix + "min_refine",0.3);
      adapt_max_coarsen[ia] = p->value (prefix + "max_coarsen",
				       0.5*adapt_min_refine[ia]);

    }

    adapt_max_level[ia] = p->value (prefix + "max_level",
				    std::numeric_limits<int>::max());

    adapt_level_exponent[ia] = p->value (prefix + "level_exponent",0.0);

    adapt_output[ia] = p->value_string (prefix + "output","");

    adapt_include_ghosts[ia] = p->value_logical (prefix + "include_ghosts",
						 false);
    const bool adapt_scheduled = 
      (p->type(prefix+"schedule:var") != parameter_unknown);

    if (adapt_scheduled) {
      p->group_set(0,"Adapt");
      p->group_push(adapt_list[ia]);
      p->group_push("schedule");
      adapt_schedule_index[ia] = read_schedule_(p, prefix);
      p->group_pop();
      p->group_pop();
    } else {
      adapt_schedule_index[ia] = -1;
    }
  }
}

//----------------------------------------------------------------------

void Config::read_balance_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Load balancing
  //--------------------------------------------------

  p->group_clear();

  balance_type = p->value_string ("Balance:type","charm");
  ASSERT1 ("Config::read_balance_",
          "Unknown Balance:type parameter %s; valid are \"charm\" or \"cello\"",
           balance_type.c_str(),
           ((balance_type == "charm") ||
            (balance_type == "cello")));

  const bool balance_scheduled = 
    (p->type("Balance:schedule:var") != parameter_unknown);

  if (balance_scheduled) {
    p->group_set(0,"Balance");
    p->group_push("schedule");
    balance_schedule_index = read_schedule_(p, "Balance");
    p->group_pop();
  } else {
    balance_schedule_index = -1;
  }
  
}  

//----------------------------------------------------------------------

void Config::read_boundary_ (Parameters * p) throw()
{

  //--------------------------------------------------
  // Boundary
  //--------------------------------------------------

  const bool multi_boundary =
    (p->type("Boundary:list") == parameter_list);
  

  num_boundary = multi_boundary ? 
    p->list_length("Boundary:list") : 1;

  boundary_list.resize(num_boundary);
  boundary_type.resize(num_boundary);
  boundary_axis.resize(num_boundary);
  boundary_face.resize(num_boundary);
  boundary_mask.resize(num_boundary);
  boundary_field_list.resize(num_boundary);

  for (int ib=0; ib<num_boundary; ib++) {

    boundary_list[ib] = multi_boundary ?
      p->list_value_string(ib,"Boundary:list","unknown") : "boundary";

    std::string prefix = "Boundary:";

    if (multi_boundary)
      prefix = prefix + boundary_list[ib] + ":";

    boundary_type[ib] = p->value_string(prefix+"type","unknown");

    std::string axis_str = p->value_string(prefix+"axis","all");
    if      (axis_str == "all") { boundary_axis[ib] = axis_all; }
    else if (axis_str == "x")   { boundary_axis[ib] = axis_x; }
    else if (axis_str == "y")   { boundary_axis[ib] = axis_y; }
    else if (axis_str == "z")   { boundary_axis[ib] = axis_z; }
    else {
      ERROR2 ("Config::read()", "Unknown %s %s",
	      (prefix+"axis").c_str(),axis_str.c_str());
    }

    std::string face_str = p->value_string(prefix+"face","all");
    if      (face_str == "all")   { boundary_face[ib] = face_all; }
    else if (face_str == "lower") { boundary_face[ib] = face_lower; }
    else if (face_str == "upper") { boundary_face[ib] = face_upper; }
    else {
      ERROR2 ("Config::read()", "Unknown %s %s",
	      (prefix+"face").c_str(),face_str.c_str());
    }

    boundary_mask[ib] = (p->type(prefix+"mask") 
			 == parameter_logical_expr);

    std::string param_str = prefix+"field_list";
    int field_list_type = p->type(param_str);
    if (field_list_type == parameter_list) {
      const int n = p->list_length(param_str);
      boundary_field_list[ib].resize(n);
      for (int index=0; index<n; index++) {
	boundary_field_list[ib][index] = p->list_value_string 
	  (index,param_str);
      }
    } else if (field_list_type == parameter_string) {
      boundary_field_list[ib].resize(1);
      boundary_field_list[ib][0] = p->value_string(param_str);
    } else if (field_list_type != parameter_unknown) {
      ERROR2 ("Config::read()", "Incorrect parameter type %d for %s",
	      field_list_type,param_str.c_str());
    }
	
  }

}

//----------------------------------------------------------------------

void Config::read_domain_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Domain
  //--------------------------------------------------

  for (int i=0; i<3; i++)  {
    domain_lower[i] = p->list_value_float(i, "Domain:lower", 0.0);
    domain_upper[i] = p->list_value_float(i, "Domain:upper", 1.0);
  }

}

//----------------------------------------------------------------------

void Config::read_field_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Field
  //--------------------------------------------------

  num_fields = p->list_length("Field:list"); 

  field_list.resize(num_fields);

  std::map<std::string, int> field_index;

  for (int i=0; i<num_fields; i++) {
    field_list[i] = p->list_value_string(i, "Field:list");
    field_index[field_list[i]] = i;
  }

  if (p->type("Field:ghost_depth") == parameter_integer) {
    field_ghost_depth[0] = p->value_integer("Field:ghost_depth",0);
    field_ghost_depth[1] = p->value_integer("Field:ghost_depth",0);
    field_ghost_depth[2] = p->value_integer("Field:ghost_depth",0);
  } else if (p->type("Field:ghost_depth") == parameter_list) {
    field_ghost_depth[0] = p->list_value_integer(0,"Field:ghost_depth",0);
    field_ghost_depth[1] = p->list_value_integer(1,"Field:ghost_depth",0);
    field_ghost_depth[2] = p->list_value_integer(2,"Field:ghost_depth",0);
  } else {
    field_ghost_depth[0] = 0;
    field_ghost_depth[1] = 0;
    field_ghost_depth[2] = 0;
  }

  field_alignment = p->value_integer("Field:alignment",8);

  field_centering[0].resize(num_fields);
  field_centering[1].resize(num_fields);
  field_centering[2].resize(num_fields);

  std::string param;

  for (int index_field=0; index_field<num_fields; index_field++) {

    param = std::string("Field:") + field_list[index_field] + ":centering";

    field_centering[0][index_field] = p->list_value_logical(0,param,true) ? 0 : 1;
    field_centering[1][index_field] = p->list_value_logical(1,param,true) ? 0 : 1;
    field_centering[2][index_field] = p->list_value_logical(2,param,true) ? 0 : 1;

  }

  // Add fields to groups (Field : <field_name> : group_list)

  field_group_list.resize(num_fields);

  for (int index_field=0; index_field<num_fields; index_field++) {

    std::string field = field_list[index_field];

    param =  std::string("Field:") + field + ":group_list";

    if (p->type(param) == parameter_list)  {
      // group_list is a list
      const int n = p->list_length(param);
      for (int i=0; i<n; i++) {
	std::string group = p->list_value_string(i,param);
	field_group_list[index_field].push_back(group);
      }
    } else if (p->type(param) == parameter_string) {
      // group_list is a string
      std::string group = p->value_string(param);
      field_group_list[index_field].push_back(group);
    }
  }

  // Add fields to groups (Group : <group_name> : field_list)

  int num_groups = p->list_length("Group:list");

  auto find_field_index = [&field_index](const std::string& field,
                                         const std::string& group) -> int
    {
      // don't use field_index[field] to access the index of a field because
      // we're not certain field_index already contains an entry for field. In
      // the case where field_index doesn't already contain an entry for field,
      // then an entry is created (with an incorrect index)

      auto search = field_index.find(field);
      if (search == field_index.end()){
        ERROR2("Config::read_field_",
               ("Can't add the \"%s\" field to the \"%s\" group because the "
                "field has not been defined"),
               field.c_str(), group.c_str());
      }
      return search->second;
    };

  for (int index_group = 0; index_group < num_groups; index_group++) {

    std::string group = p->list_value_string(index_group, "Group:list") ;

    param = std::string("Group:") + group + ":field_list";

    if (p->type(param) == parameter_list) {
      // field_list is a list
      const int n = p->list_length(param);
      for (int i=0; i<n; i++) {
	std::string field = p->list_value_string(i,param);
	const int index_field = find_field_index(field,group);
	field_group_list[index_field].push_back(group);
      }
    } else if (p->type(param) == parameter_string) {
      // field_list is a string
      std::string field = p->value_string(param);
      const int index_field = find_field_index(field,group);
      field_group_list[index_field].push_back(group);
    }
  }

  field_padding = p->value_integer("Field:padding",0);

  field_history = p->value_integer("Field:history",0);

  // Field precision

  std::string precision_str = p->value_string("Field:precision","default");

  if      (precision_str == "default")   field_precision = precision_default;
  else if (precision_str == "single")    field_precision = precision_single;
  else if (precision_str == "double")    field_precision = precision_double;
  else if (precision_str == "quadruple") field_precision = precision_quadruple;
  else {
    ERROR1 ("Config::read()", "Unknown precision %s",
	    precision_str.c_str());
  }

  field_prolong   = p->value_string ("Field:prolong","enzo");
  field_restrict  = p->value_string ("Field:restrict","linear");
}

//----------------------------------------------------------------------

void Config::read_initial_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Initial
  //--------------------------------------------------

  TRACE("Parameters: Initial");
  initial_cycle = p->value_integer("Initial:cycle",0);
  initial_time  = p->value_float  ("Initial:time",0.0);
  initial_new = p->value_logical ("Initial:new",false);

  num_initial = p->list_length("Initial:list");

  initial_list.resize(num_initial);
  for (int i=0; i<num_initial; i++) {

    std::string name = 
      p->list_value_string(i,"Initial:list");

    initial_list[i] = name;

  }

  // Restart

  initial_restart      = p->value_logical ("Initial:restart",false);
  initial_restart_dir  = p->value_string  ("Initial:restart_dir","");

  // InitialTrace
  initial_trace_name = p->value_string ("Initial:trace:name","trace");
  initial_trace_field = p->value_string ("Initial:trace:field","");
  initial_trace_mpp = p->value_float ("Initial:trace:mass_per_particle",0.0);
  initial_trace_dx = p->list_value_integer (0,"Initial:trace:stride",1);
  initial_trace_dy = p->list_value_integer (1,"Initial:trace:stride",1);
  initial_trace_dz = p->list_value_integer (2,"Initial:trace:stride",1);
}

//----------------------------------------------------------------------

void Config::read_memory_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Memory
  //--------------------------------------------------

  memory_active = p->value_logical("Memory:active",true);
  memory_warning_mb =  p->value_float("Memory:warning_mb",0.0);
  memory_limit_gb =    p->value_float("Memory:limit_gb",0.0);
}

//----------------------------------------------------------------------

void Config::read_mesh_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Mesh
  //--------------------------------------------------

  TRACE("Parameters: Mesh");
  mesh_root_rank = p->value_integer("Mesh:root_rank",0);

  // Adjust ghost zones for unused ranks
  if (mesh_root_rank < 2) field_ghost_depth[1] = 0;
  if (mesh_root_rank < 3) field_ghost_depth[2] = 0;
  
  //--------------------------------------------------

  int mx = mesh_root_blocks[0] = p->list_value_integer(0,"Mesh:root_blocks",1);
  int my = mesh_root_blocks[1] = p->list_value_integer(1,"Mesh:root_blocks",1);
  int mz = mesh_root_blocks[2] = p->list_value_integer(2,"Mesh:root_blocks",1);

  const int m = mx*my*mz;

  if ( ! (m >= CkNumPes()) ) {
    WARNING4 ("Config::read_mesh_()",
	    "Number of root blocks %d x %d x %d cannot be be "
	    "less than number of processes %d",
	    mx,my,mz,CkNumPes());
  }

  //--------------------------------------------------

  mesh_root_size[0] = p->list_value_integer(0,"Mesh:root_size",1);
  mesh_root_size[1] = p->list_value_integer(1,"Mesh:root_size",1);
  mesh_root_size[2] = p->list_value_integer(2,"Mesh:root_size",1);

  //--------------------------------------------------

  mesh_max_level = p->value_integer
    ("Adapt:max_level",0);
  mesh_max_initial_level = p->value_integer
    ("Adapt:max_initial_level",mesh_max_level);

  // Note mesh_min_level may be < 0 for multigrid

  mesh_min_level = p->value_integer("Adapt:min_level",0);

  if ( mesh_min_level > 0 ) {
    ERROR1 ("Config::read", 
		    "The value of mesh_min_level: %d should be less than or equal to zero", 
		    mesh_min_level);
  }

  // Handle 1D and 2D simulations by adjusting the number of cells along the extra dimensions
  if (mesh_root_rank < 2) mesh_root_size[1] = 1;
  if (mesh_root_rank < 3) mesh_root_size[2] = 1;

  // Dimensions of the active zone on each block along each axis
  int ax = mesh_root_size[0] / mesh_root_blocks[0]; 
  int ay = mesh_root_size[1] / mesh_root_blocks[1];
  int az = mesh_root_size[2] / mesh_root_blocks[2];

  //  Constraints on the block size based on the ghost depth
  if ( mesh_max_level > 0 ) {
    if ( !(ax >= 2*field_ghost_depth[0] && ay >= 2*field_ghost_depth[1] && az >= 2*field_ghost_depth[2] ) ) {
      ERROR3 ("Config::read", 
		"Dimensions of the active zone on each block (%d, %d, %d) should be at least double the size of the ghost depth for AMR simulations: ", 
		ax, ay, az);
    }  
    if ( (ax%2 != 0) || (ay%2 != 0) && (az%2 != 0) ) {
      ERROR3 ("Config::read",
  		      "Dimensions of the active zone on each block (%d, %d, %d) should each be even for AMR simulations" ,
		      ax, ay, az);
    }  
  }
  else if ( mesh_max_level == 0 ) {   
    if ( !(ax >= field_ghost_depth[0] && ay >= field_ghost_depth[1] && az >= field_ghost_depth[2] ) ) {
      ERROR3 ("Config::read",
  		      "Dimensions of the active zone on each block (%d, %d, %d) should be at least as large as the ghost depth",
		      ax, ay, az);
    }  
  }
}

//----------------------------------------------------------------------

void Config::read_method_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Method
  //--------------------------------------------------

  TRACE("Parameters: Method");

  num_method = p->list_length("Method:list");

  method_list.   resize(num_method);
  method_courant.resize(num_method);
  method_file_name.resize(num_method);
  method_path_name.resize(num_method);
  method_debug_print.resize(num_method);
  method_debug_coarse.resize(num_method);
  method_debug_ghost.resize(num_method);
  method_flux_correct_group.resize(num_method);
  method_flux_correct_enable.resize(num_method);
  method_flux_correct_min_digits_fields.resize(num_method);
  method_flux_correct_min_digits_values.resize(num_method);
  method_field_list.resize(num_method);
  method_particle_list.resize(num_method);
  method_output_blocking[0].resize(num_method);
  method_output_blocking[1].resize(num_method);
  method_output_blocking[2].resize(num_method);
  method_output_all_blocks.resize(num_method);
  method_prolong.resize(num_method);
  method_ghost_depth.resize(num_method);
  method_min_face_rank.resize(num_method);
  method_all_fields.resize(num_method);
  method_all_particles.resize(num_method);
  method_timestep.resize(num_method);
  method_schedule_index.resize(num_method);
  method_close_files_seconds_stagger.resize(num_method);
  method_close_files_seconds_delay.resize(num_method);
  method_close_files_group_size.resize(num_method);
  method_trace_name.resize(num_method);
  method_type.resize(num_method);
  
  method_courant_global = p->value_float ("Method:courant",1.0);
  
  for (int index_method=0; index_method<num_method; index_method++) {

    std::string name = 
      p->list_value_string(index_method,"Method:list");

    std::string full_name = std::string("Method:") + name;

    method_list[index_method] = name;

    // Read schedule for the Method object if any
      
    const bool method_scheduled = 
      (p->type(full_name + ":schedule:var") != parameter_unknown);

    if (method_scheduled) {
      p->group_set(0,"Method");
      p->group_push(name);
      p->group_push("schedule");
      method_schedule_index[index_method] = read_schedule_(p, name);
      p->group_pop();
    } else {
      method_schedule_index[index_method] = -1;
    }

    // Read method file_name
    if (p->type(full_name+":file_name") == parameter_string) {
      method_file_name[index_method].push_back
        (p->value_string(full_name+":file_name",""));
    } else if (p->type(full_name+":file_name") == parameter_list) {
      int size = p->list_length(full_name+":file_name");
      if (size > 0) method_file_name[index_method].resize(size);
      for (int i=0; i<size; i++) {
        method_file_name[index_method][i] =
          p->list_value_string(i,full_name+":file_name","");
      }
    }

    // Read method path_name
    if (p->type(full_name+":path_name") == parameter_string) {
      method_path_name[index_method].push_back
        (p->value_string(full_name+":path_name",""));
    } else if (p->type(full_name+":path_name") == parameter_list) {
      int size = p->list_length(full_name+":path_name");
      if (size > 0) method_path_name[index_method].resize(size);
      for (int i=0; i<size; i++) {
        method_path_name[index_method][i] =
          p->list_value_string(i,full_name+":path_name","");
      }
    }

    // Read throttling parameters for MethodCloseFiles
    method_close_files_seconds_stagger[index_method] = p->value_float
      (full_name + ":seconds_stagger",0.0);
    method_close_files_seconds_delay[index_method] = p->value_float
      (full_name + ":seconds_delay",0.0);
    method_close_files_group_size[index_method] = p->value_integer
      (full_name + ":group_size",std::numeric_limits<int>::max());

    // Read courant condition if any
    method_courant[index_method] = p->value_float  (full_name + ":courant",1.0);

    // Read any MethodDebug parameters
    method_debug_print[index_method] = p->value_logical
      (full_name + ":print",false);
    method_debug_coarse[index_method] = p->value_logical
      (full_name + ":coarse",false);
    method_debug_ghost[index_method] = p->value_logical
      (full_name + ":ghost",false);

    // Read field group for flux correction
    method_flux_correct_group[index_method] =
      p->value_string (full_name + ":group","conserved");
    method_flux_correct_enable[index_method] =
      p->value_logical (full_name + ":enable",true);

    std::string min_digits_name = full_name + ":min_digits";
    if (p->type(min_digits_name) == parameter_float){
      // backwards compatibility
      method_flux_correct_min_digits_fields[index_method] = {"density"};
      method_flux_correct_min_digits_values[index_method].push_back
        (p->value_float (min_digits_name, 0.0));
    } else if (p->type(min_digits_name) == parameter_list){
      // load pairs of fields and min_digits
      int list_length = p->list_length(min_digits_name);
      ASSERT1("Config::read",
              "The list assigned to %s must have a non-negative, even length",
              min_digits_name.c_str(),
              (list_length >= 0) && (list_length % 2 == 0));
      for (int i =0; i < list_length; i+=2){
        method_flux_correct_min_digits_fields[index_method].push_back
          (p->list_value_string(i, min_digits_name));
        method_flux_correct_min_digits_values[index_method].push_back
          (p->list_value_float (i+1, min_digits_name, 0.0));
      }
    } else if (p->param(min_digits_name) != nullptr){
      ERROR1("Config::read", "%s has an invalid type", min_digits_name.c_str());
    }

    method_flux_correct_single_array =
      p->value_logical (full_name + ":single_array",true);

    // Field and particle lists if needed by MethodRefresh
    int n = p->list_length(full_name + ":field_list");
    method_field_list[index_method].resize(n);
    for (int i=0; i<n; i++) {
      method_field_list[index_method][i] =
        p->list_value_string(i,full_name+":field_list");
    }
    n = p->list_length(full_name + ":particle_list");
    method_particle_list[index_method].resize(n);
    for (int i=0; i<n; i++) {
      method_particle_list[index_method][i] =
        p->list_value_string(i,full_name+":particle_list");
    }

    for (int i=0; i<3; i++) {
      method_output_blocking[i][index_method] =
        p->list_value_integer(i,full_name+":blocking",1);
    }
    method_output_all_blocks[index_method] =
      p->value_logical(full_name+":all_blocks",true);

    method_prolong[index_method] =
      p->value_string(full_name+":prolong","linear");

    // Read refresh method parameters
    method_ghost_depth[index_method] =
      p->value_integer(full_name+":ghost_depth",0);
    method_min_face_rank[index_method] =
      p->value_integer(full_name+":min_face_rank",0); // default 0 all faces
    method_all_fields[index_method] =
      p->value_logical(full_name+":all_fields",false);
    method_all_particles[index_method] =
      p->value_logical(full_name+":all_particles",false);

    // Read specified timestep, if any (for MethodTrace)
    method_timestep[index_method] = p->value_float
      (full_name + ":timestep",std::numeric_limits<double>::max());

    method_trace_name[index_method] = p->value_string
      (full_name + ":name", "trace");

    method_type[index_method] = p->value_string
      (full_name + ":type", name);
  }
  method_null_dt = p->value_float
    ("Method:null:dt",std::numeric_limits<double>::max());

}

//----------------------------------------------------------------------

void Config::read_monitor_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Monitor
  //--------------------------------------------------

  monitor_debug   = p->value_logical("Monitor:debug",  false);
  monitor_verbose = p->value_logical("Monitor:verbose",false);

}

//----------------------------------------------------------------------

void Config::read_output_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Output
  //--------------------------------------------------

  p->group_set(0,"Output");

  num_output = p->list_length("list");

  output_list.resize(num_output);
  output_type.resize(num_output);
  output_axis.resize(num_output);
  output_colormap.resize(num_output);
  output_image_lower.resize(num_output);
  output_image_upper.resize(num_output);
  output_image_type.resize(num_output);
  output_image_log.resize(num_output);
  output_image_abs.resize(num_output);
  output_image_mesh_color.resize(num_output);
  output_image_mesh_order.resize(num_output);
  output_image_color_particle_attribute.resize(num_output);
  output_image_size.resize(num_output);
  output_image_reduce_type.resize(num_output);
  output_image_ghost.resize(num_output);
  output_image_face_rank.resize(num_output);
  output_image_min.resize(num_output);
  output_image_max.resize(num_output);
  output_min_level.resize(num_output);
  output_max_level.resize(num_output);
  output_leaf_only.resize(num_output);
  output_schedule_index.resize(num_output);
  output_dir.resize(num_output);
  output_stride_write.resize(num_output);
  output_stride_wait.resize(num_output);
  output_field_list.resize(num_output);
  output_particle_list.resize(num_output);
  output_name.resize(num_output);

  output_dir_global = p->value_string("dir_global",".");

  for (int index_output=0; index_output<num_output; index_output++) {

    TRACE1 ("index = %d",index_output);

    output_list[index_output] = 
      p->list_value_string (index_output,"Output:list","unknown");

    p->group_set(1,output_list[index_output]);

    output_type[index_output] = p->value_string("type","unknown");

    if (output_type[index_output] == "unknown") {
      ERROR1("Config::read",
	     "Output:%s:type parameter is undefined",
	     output_list[index_output].c_str());
    }

    output_stride_write[index_output] = p->value_integer("stride_write",0);

    output_stride_wait[index_output] = p->value_integer("stride_wait",0);

    if (p->type("dir") == parameter_string) {
      output_dir[index_output].resize(1);
      output_dir[index_output][0] = p->value_string("dir","");
    } else if (p->type("dir") == parameter_list) {
      int size = p->list_length("dir");
      if (size > 0) output_dir[index_output].resize(size);
      for (int i=0; i<size; i++) {
	output_dir[index_output][i] = p->list_value_string(i,"dir","");
      }
    }

    if (p->type("name") == parameter_string) {
      output_name[index_output].resize(1);
      output_name[index_output][0] = p->value_string("name","");
    } else if (p->type("name") == parameter_list) {
      int size = p->list_length("name");
      if (size > 0) output_name[index_output].resize(size);
      for (int i=0; i<size; i++) {
	output_name[index_output][i] = p->list_value_string(i,"name","");
      }
    }

    // // File group (data dump only)
    // if (p->type("group") == parameter_string) {
    //   output_group[index_output].resize(1);
    //   output_group[index_output][0] = p->value_string("group","");
    // } else if (p->type("group") == parameter_list) {
    //   int size = p->list_length("group");
    //   if (size > 0) output_group[index_output].resize(size);
    //   for (int i=0; i<size; i++) {
    // 	output_group[index_output][i] = p->list_value_string(i,"group","");
    //   }
    // }

    if (p->type("field_list") == parameter_list) {
      int length = p->list_length("field_list");
      output_field_list[index_output].resize(length);
      for (int i=0; i<length; i++) {
	output_field_list[index_output][i] = p->list_value_string(i,"field_list","");
      }
    }

    if (p->type("particle_list") == parameter_list) {
      int length = p->list_length("particle_list");
      output_particle_list[index_output].resize(length);
      for (int i=0; i<length; i++) {
	output_particle_list[index_output][i] = p->list_value_string(i,"particle_list","");
      }
    }

    // Read schedule for the Output object
      
    p->group_push("schedule");
    output_schedule_index[index_output] = 
      read_schedule_(p, output_list[index_output]);
    p->group_pop();

    // Image 
    
    if (output_type[index_output] == "image") {


      if (p->type("axis") != parameter_unknown) {
	std::string axis = p->value_string("axis");
	ASSERT2("Problem::initialize_output",
		"Output %s axis %s must be \"x\", \"y\", or \"z\"",
		output_list[index_output].c_str(), axis.c_str(),
		axis=="x" || axis=="y" || axis=="z");
	output_axis[index_output] = axis;
      }  else {
	WARNING1 ("Config::read()",
		  "output_axis[%d] set to z",index_output);

	output_axis[index_output] = "z";
      }


      output_image_type[index_output] = p->value_string("image_type","data");

      output_image_log[index_output] = p->value_logical("image_log",false);
      output_image_abs[index_output] = p->value_logical("image_abs",false);

      output_image_mesh_color[index_output] = 
	p->value_string("image_mesh_color","level");
      output_image_mesh_order[index_output] = 
	p->value_string("image_mesh_order","none");

      output_image_color_particle_attribute[index_output] = 
	p->value_string("image_color_particle_attribute","");

      output_image_size[index_output].resize(2);
      output_image_size[index_output][0] = 
	p->list_value_integer(0,"image_size",512);
      output_image_size[index_output][1] = 
	p->list_value_integer(1,"image_size",512);

      output_image_reduce_type[index_output] = 
	p->value_string("image_reduce_type","sum");

      output_image_face_rank[index_output] = 
	p->value_integer("image_face_rank",3);

      output_image_ghost[index_output] = 
	p->value_logical("image_ghost",false);

      output_image_min[index_output] =
	p->value_float("image_min",std::numeric_limits<double>::max());
      output_image_max[index_output] =
	p->value_float("image_max",-std::numeric_limits<double>::max());

      output_min_level[index_output] = p->value_integer("min_level",0);
      output_max_level[index_output] =
	p->value_integer("max_level",std::numeric_limits<int>::max());
      output_leaf_only[index_output] = p->value_logical("leaf_only",true);

      if (p->type("colormap") == parameter_list) {
	int size = p->list_length("colormap");
	output_colormap[index_output].resize(size);
	for (int i=0; i<size; i++) {
	  output_colormap[index_output][i] = 
	    p->list_value_float(i,"colormap",0.0);
	}
      }

      output_image_lower[index_output].resize(3);
      output_image_upper[index_output].resize(3);
      for (int axis=0; axis<3; axis++) {
	output_image_lower[index_output][axis] =
	  p->list_value_float(axis,"image_lower",-std::numeric_limits<double>::max());
	output_image_upper[index_output][axis] =
	  p->list_value_float(axis,"image_upper",std::numeric_limits<double>::max());
      }

    }
  }  

}

//----------------------------------------------------------------------

void Config::read_particle_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Particle
  //--------------------------------------------------

  particle_batch_size = p->value_integer("Particle:batch_size",1024);

  num_particles = p->list_length("Particle:list"); 

  particle_list.resize(num_particles);
  particle_interleaved.resize(num_particles);
  particle_constant_name.resize(num_particles);
  particle_constant_type.resize(num_particles);
  particle_constant_value.resize(num_particles);
  particle_attribute_name.resize(num_particles);
  particle_attribute_type.resize(num_particles);
  particle_attribute_position[0].resize(num_particles);
  particle_attribute_position[1].resize(num_particles);
  particle_attribute_position[2].resize(num_particles);
  particle_attribute_velocity[0].resize(num_particles);
  particle_attribute_velocity[1].resize(num_particles);
  particle_attribute_velocity[2].resize(num_particles);

  // ... first map attribute scalar type name to type_enum int
  std::map<std::string,int> type_val;
  for (int i=0; i<NUM_TYPES; i++) {
    type_val[cello::type_name[i]] = i;
  }

  for (int it=0; it<num_particles; it++) {

    std::string name_type = p->list_value_string(it, "Particle:list");

    particle_list [it]   = name_type;
    particle_index[name_type] = it;

    std::string type_str = "Particle:" + name_type;

    // are attributes are interleaved?

    particle_interleaved[it] = 
      p->value_logical(type_str+":interleaved",false);

    // Particle:<type>:constants list elements contain name, type, and
    // value

    std::string const_str = type_str + ":constants";
    
    const int nc3 = p->list_length(const_str);

    ASSERT2 ("read_particle_",
	     "Particle type %d constants list length %d "
	     "must be divisible by three",
	     it,nc3,nc3%3==0);

    const int nc = nc3 / 3;

    particle_constant_name[it].resize(nc);
    particle_constant_type[it].resize(nc);
    particle_constant_value[it].resize(nc);

    for (int ia=0; ia<nc; ia++) {

      std::string name = p->list_value_string (3*ia,  const_str,"unknown");
      std::string type = p->list_value_string (3*ia+1,const_str,"unknown");

      if (type == "default") {
        type = default_precision_string;
      }

      ASSERT3 ("read_particle_",
	       "Particle type %d constant %d has unknown constant name %s",
	       it,ia,name.c_str(),
	       name != "unknown");

      particle_constant_name[it][ia]  = name;
      particle_constant_type[it][ia]  = type;

      if (cello::type_is_float(type_val[type])) {
	particle_constant_value[it][ia] = 
	  p->list_value_float (3*ia+2,const_str,0.0);
      } else if (cello::type_is_int(type_val[type])) {
	particle_constant_value[it][ia] =
	  p->list_value_integer (3*ia+2,const_str,0);
      }

    }

    // Particle:<type>:attributes list elements alternate attribute
    // name and its type (see type_enum in cello.hpp)

    std::string attrib_str = type_str + ":attributes";
    
    const int na2 = p->list_length(attrib_str);

    ASSERT1 ("read_particle_",
	     "Particle type %d has no attributes",
	     it,na2>0);

    ASSERT2 ("read_particle_",
	     "Particle type %d attributes list length %d must be even",
	     it,na2,na2%2==0);

    const int na = na2 / 2;

    particle_attribute_name[it].resize(na);
    particle_attribute_type[it].resize(na);

    std::map <std::string,int> attribute_index;

    for (int ia=0; ia<na; ia++) {

      std::string name = p->list_value_string (2*ia  ,attrib_str,"unknown");
      std::string type = p->list_value_string (2*ia+1,attrib_str,"unknown");

      ASSERT3 ("read_particle_",
	       "Particle type %d attribute %d has unknown attribute name %s",
	       it,ia,name.c_str(),
	       name != "unknown");

      particle_attribute_name[it][ia]  = name;
      particle_attribute_type[it][ia]  = type;
     
      attribute_index[name] = ia;
    }

    // Particle:<type>:position
    // Particle:<type>:velocity

    const std::string param_position = type_str+":position";
    const std::string param_velocity = type_str+":velocity";

    for (int axis=0; axis<3; axis++) {

      std::string pos = p->list_value_string (axis,param_position,"");
      int ia_p = (pos != "") ? attribute_index[pos] : -1;
      particle_attribute_position[axis][it] = ia_p;

      std::string vel = p->list_value_string (axis,param_velocity,"");
      int ia_v = (vel != "") ? attribute_index[vel] : -1;
      particle_attribute_velocity[axis][it] = ia_v;
    }
  }

  // Add particles to groups (Particle : <particle_name> : group_list)

  particle_group_list.resize(num_particles);

  for (int index_particle=0; index_particle<num_particles; index_particle++) {

    std::string particle = particle_list[index_particle];

    std::string param =  std::string("Particle:") + particle + ":group_list";

    if (p->type(param) == parameter_list)  {
      // group_list is a list
      const int n = p->list_length(param);
      for (int i=0; i<n; i++) {
	std::string group = p->list_value_string(i,param);
	particle_group_list[index_particle].push_back(group);
      }
    } else if (p->type(param) == parameter_string) {
      // group_list is a string
      std::string group = p->value_string(param);
      particle_group_list[index_particle].push_back(group);
    }
  }

  // Add particles to groups (Group : <group_name> : particle_list)

  int num_groups = p->list_length("Group:list"); 

  for (int index_group = 0; index_group < num_groups; index_group++) {

    std::string group = p->list_value_string(index_group, "Group:list") ;

    std::string param = std::string("Group:") + group + ":particle_list";

    if (p->type(param) == parameter_list) {
      // particle_list is a list
      const int n = p->list_length(param);
      for (int i=0; i<n; i++) {
	std::string particle = p->list_value_string(i,param);
	const int index_particle = particle_index[particle];
	particle_group_list[index_particle].push_back(group);
      }
    } else if (p->type(param) == parameter_string) {
      // particle_list is a string
      std::string particle = p->value_string(param);
      const int index_particle = particle_index[particle];
      particle_group_list[index_particle].push_back(group);
    }
  }
}

//----------------------------------------------------------------------

void Config::read_performance_ (Parameters * p) throw()
{
#ifdef CONFIG_USE_PAPI  
  if (p->type("Performance:papi:counters") == parameter_list) {
    int length = p->list_length("Performance:papi:counters");
    performance_papi_counters.resize(length);
    for (int i=0; i<length; i++) {
      performance_papi_counters[i] = p->list_value_string
	(i,"Performance:papi:counters","");
      TRACE2("performance_papi_counters[%d] = %s",
	     i,performance_papi_counters[i].c_str());
    }
  }
#endif  

  performance_warnings = p->value_logical("Performance:warnings",false);

#ifdef CONFIG_USE_PROJECTIONS
  
  int i_on = -1;
  int i_off = -1;
  
  if (p->type("Performance:projections:schedule_on:var") != parameter_unknown) {
    p->group_set(0,"Performance");
    p->group_push("projections");
    p->group_push("schedule_on");
    i_on = read_schedule_(p,"projections_on");
  }
  if (p->type("Performance:projections:schedule_off:var") != parameter_unknown) {
    p->group_set(0,"Performance");
    p->group_push("projections");
    p->group_push("schedule_off");
    i_off = read_schedule_(p,"projections_off");
  }
  p->group_clear();

  performance_projections_on_at_start =  p->value_logical
    ("Performance:projections:on_at_start",true);
  
  // Check that both projections_on and off schedules are defined or undefined together
  if ((i_on == -1 && i_off == -1) || (i_on != -1 && i_off != -1)) {
    performance_on_schedule_index  = i_on;
    performance_off_schedule_index = i_off;
  } else {
    ERROR2("Config::read_performance-()",
	   "Performance:projections:schedule_on [%d] and Performance:projections:schedule_off [%d]\n"
	   "must be both defined or both undefined",
	   i_on,i_off);
  }
#endif    

}

//----------------------------------------------------------------------

void Config::read_physics_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Physics
  //--------------------------------------------------

  num_physics = p->list_length("Physics:list"); 

  physics_list.resize(num_physics);
  
  for (int index_physics=0; index_physics<num_physics; index_physics++) {

    std::string name = 
      p->list_value_string(index_physics,"Physics:list");

    physics_list[index_physics] = name;

  }
}

//----------------------------------------------------------------------

void Config::read_solver_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Solver
  //--------------------------------------------------

  TRACE("Parameters: Solver");

  num_solvers = p->list_length("Solver:list");

  solver_list         .resize(num_solvers);
  solver_type         .resize(num_solvers);
  solver_solve_type   .resize(num_solvers);
  solver_iter_max     .resize(num_solvers);
  solver_res_tol      .resize(num_solvers);
  solver_diag_precon  .resize(num_solvers);
  solver_monitor_iter .resize(num_solvers);
  solver_restrict     .resize(num_solvers);
  solver_prolong      .resize(num_solvers);
  solver_min_level    .resize(num_solvers);
  solver_max_level    .resize(num_solvers);
  solver_field_x      .resize(num_solvers);
  solver_field_b      .resize(num_solvers);

  for (int index_solver=0; index_solver<num_solvers; index_solver++) {

    std::string name = 
      p->list_value_string(index_solver,"Solver:list");

    std::string full_name = std::string("Solver:") + name;

    solver_list[index_solver] = name;

    solver_index[name] = index_solver;

    solver_type[index_solver] = p->value_string (full_name + ":type","unknown");
    
    solver_type[index_solver] = p->value_string (full_name + ":type","unknown");

    solver_solve_type[index_solver] = p->value_string
      (full_name + ":solve_type","leaf");

    solver_iter_max[index_solver] = p->value_integer
      (full_name + ":iter_max",1000);
    
    solver_res_tol[index_solver] = p->value_float
      (full_name + ":res_tol",1e-6);

    solver_diag_precon[index_solver] = p->value_logical
      (full_name + ":diag_precon",false);
    
    solver_monitor_iter[index_solver] = p->value_integer
      (full_name + ":monitor_iter",0);

    solver_restrict[index_solver] = p->value_string
      (full_name + ":restrict","linear");

    solver_prolong[index_solver] = p->value_string
      (full_name + ":prolong","linear");

    solver_min_level[index_solver] = p->value_integer
      (full_name + ":min_level",0);

    solver_max_level[index_solver] = p->value_integer
      (full_name + ":max_level",mesh_max_level);

    solver_field_x[index_solver] = p->value_string
      (full_name + ":field_x","unknown");

    solver_field_b[index_solver] = p->value_string
      (full_name + ":field_b","unknown");
  }  
}

//----------------------------------------------------------------------

void Config::read_stopping_ (Parameters * p) throw()
{
  stopping_cycle = p->value_integer
    ( "Stopping:cycle" , std::numeric_limits<int>::max() );
  stopping_time  = p->value_float
    ( "Stopping:time" , std::numeric_limits<double>::max() );

  stopping_seconds = std::numeric_limits<double>::max();

  if (p->type("Stopping:seconds") != parameter_unknown) {
    stopping_seconds  = p->value_float
      ( "Stopping:seconds" , std::numeric_limits<double>::max() );
  } else if (p->type("Stopping:minutes") != parameter_unknown) {
    stopping_seconds  = p->value_float
      ( "Stopping:minutes" , std::numeric_limits<double>::max() );
    stopping_seconds *= 60;
  } else if (p->type("Stopping:hours") != parameter_unknown) {
    stopping_seconds  = p->value_float
      ( "Stopping:hours" , std::numeric_limits<double>::max() );
    stopping_seconds *= 3600;
  }

  stopping_interval = p->value_integer ( "Stopping:interval" , 1);
}

void Config::read_units_ (Parameters * p) throw()
{
  //======================================================================
  // Units
  //======================================================================
  
  units_mass    = p->value_float ("Units:mass",   1.0);
  units_density = p->value_float ("Units:density",1.0);
  units_length  = p->value_float ("Units:length", 1.0);
  units_time    = p->value_float ("Units:time",   1.0);

}

//----------------------------------------------------------------------

void Config::read_testing_ (Parameters * p) throw()
{
  testing_cycle_final = p->value_integer("Testing:cycle_final",0);
  if (p->type("Testing:time_final") == parameter_list) {
    int length = p->list_length("Testing:time_final");
    testing_time_final.resize(length);
    for (int i=0; i<length; i++) {
      testing_time_final[i] = p->list_value_float (i,"Testing:time_final",0.0);
    }
  } else {
    testing_time_final.resize(1);
    testing_time_final[0]  = p->value_float  ("Testing:time_final", 0.0);
  }
  testing_time_tolerance = p->value_float  ("Testing:time_tolerance", 1e-6);
}

//======================================================================

int Config::read_schedule_(Parameters * p, const std::string group)
{
  int index = index_schedule;

  schedule_list.resize(index+1);
  schedule_type.resize(index+1);
  schedule_var.resize(index+1);
  schedule_start.resize(index+1);
  schedule_stop.resize(index+1);
  schedule_step.resize(index+1);
  
  std::string var = p->value_string("var","none");

  schedule_var[index] = var;

  bool var_is_int = true;

  // Get variable associated with the schedule 
  if      (schedule_var[index] == "cycle")    var_is_int = true;
  else if (schedule_var[index] == "time")     var_is_int = false;
  else if (schedule_var[index] == "seconds")  var_is_int = false;
  else if (schedule_var[index] == "minutes")  var_is_int = false;
  else if (schedule_var[index] == "hours")    var_is_int = false;
  else {
    ERROR2 ("Config::read",
	    "Schedule variable %s is not recognized for parameter group %s",
	    schedule_var[index].c_str(),group.c_str());
  }

  // Determine the schedule type (interval or list)

  const bool type_is_list = (p->type("list") != parameter_unknown);
  const bool type_is_interval = ! type_is_list;

  if (type_is_interval) schedule_type[index] = "interval";
  if (type_is_list)     schedule_type[index] = "list";

  const int    max_int    = std::numeric_limits<int>::max();
  const double max_double = std::numeric_limits<double>::max();

  if (type_is_interval) {
    if (var_is_int) {
      schedule_start[index] = p->value("start",0);
      schedule_step[index]  = p->value("step",1);
      schedule_stop[index]  = p->value("stop",max_int);
    } else {
      schedule_start[index] = p->value("start",0.0);
      schedule_step[index]  = p->value("step",1.0);
      schedule_stop[index]  = p->value("stop",max_double);
    }
  } else if (type_is_list) {
    int n = p->list_length("list");
    if (n == 0) {
      ERROR1 ("Config::read",
	      "Schedule variable %s has length 0",
	      (group + ":list").c_str());
    }
    schedule_list[index].resize(n);
    if (var_is_int) {
      for (int i=0; i<n; i++) {
	schedule_list[index][i] = p->value(i,"list",0);
      }
    } else {
      for (int i=0; i<n; i++) {
	schedule_list[index][i] = p->value(i,"list",0.0);
      }
    }
  } else {
    ERROR2 ("Config::read_schedule_",
	    "Schedule type %s is not recognized for parameter group %s",
	    schedule_type[index].c_str(),group.c_str());
  }

  return index_schedule++;
}
//======================================================================

