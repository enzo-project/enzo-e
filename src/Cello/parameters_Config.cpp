// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Config.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-03
/// @brief    Implementation of the Config class 

#include "cello.hpp"
#include "parameters.hpp"

//----------------------------------------------------------------------

void Config::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
  p | num_boundary;
  PUParray(p,boundary_list,MAX_BOUNDARY);
  PUParray(p,boundary_type,MAX_BOUNDARY);
  PUParray(p,boundary_axis,MAX_BOUNDARY);
  PUParray(p,boundary_face,MAX_BOUNDARY);
  PUParray(p,boundary_mask,MAX_BOUNDARY);
  PUParray(p,boundary_field_list,MAX_BOUNDARY);

  PUParray(p,domain_lower,3);
  PUParray(p,domain_upper,3);

  p | num_fields;
  p | field_alignment;
  PUParray(p,field_centering,3);
  p | field_courant;
  p | field_fields;
  PUParray(p,field_ghosts,3);
  if (p.isUnpacking()) {
    TRACE3 ("gx,gy,gz = %d %d %d",field_ghosts[0],field_ghosts[1],field_ghosts[2]);
  }

  p | field_padding;
  p | field_precision;
  p | field_refresh_rank;

  p | initial_cycle;
  p | initial_type;
  p | initial_time;
  p | initial_max_level;

  p | memory_active;

  PUParray(p,mesh_root_blocks,3);
  p | mesh_root_rank;
  PUParray(p,mesh_root_size,3);
  p | mesh_max_level;
  p | mesh_balance;
  p | mesh_adapt_type;
  p | mesh_adapt_fields;
  p | mesh_adapt_slope_min_refine;
  p | mesh_adapt_slope_max_coarsen;
  p | mesh_adapt_mass_min;
  p | mesh_adapt_mass_level_exponent;
  p | mesh_adapt_mass_min_overdensity;
  p | mesh_adapt_balance;
  p | mesh_adapt_interval;

  p | control_sync_adapt_enter;
  p | control_sync_adapt_called;
  p | control_sync_adapt_end;
  p | control_sync_adapt_exit;
  p | control_sync_adapt_next;
  p | control_sync_compute_enter;
  p | control_sync_compute_exit;
  p | control_sync_exit;
  p | control_sync_output_enter;
  p | control_sync_output_exit;
  p | control_sync_refresh_enter;
  p | control_sync_refresh_exit;
  p | control_sync_stopping_enter;
  p | control_sync_stopping_exit;

  p | method_sequence;

  p | monitor_debug;

  p | num_file_groups;
  p | output_file_groups;
  PUParray (p,output_type,MAX_FILE_GROUPS);
  PUParray (p,output_image_axis,MAX_FILE_GROUPS);
  PUParray (p,output_image_block_size,MAX_FILE_GROUPS);
  PUParray (p,output_image_colormap,MAX_FILE_GROUPS);
  PUParray (p,output_image_type,MAX_FILE_GROUPS);
  PUParray (p,output_image_log,MAX_FILE_GROUPS);
  PUParray (p,output_image_mesh_color,MAX_FILE_GROUPS);
  PUParray (p,output_image_size,MAX_FILE_GROUPS);
  PUParray (p,output_image_reduce_type,MAX_FILE_GROUPS);
  PUParray (p,output_image_ghost,MAX_FILE_GROUPS);
  PUParray (p,output_image_specify_bounds,MAX_FILE_GROUPS);
  PUParray (p,output_image_min,MAX_FILE_GROUPS);
  PUParray (p,output_image_max,MAX_FILE_GROUPS);
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

  p | performance_papi_counters;
  p | performance_name;
  p | performance_stride;

  p | projections_schedule_on_type;
  p | projections_schedule_on_var;
  p | projections_schedule_on_start;
  p | projections_schedule_on_stop;
  p | projections_schedule_on_step;
  p | projections_schedule_on_list;
  p | projections_schedule_off_type;
  p | projections_schedule_off_var;
  p | projections_schedule_off_start;
  p | projections_schedule_off_stop;
  p | projections_schedule_off_step;
  p | projections_schedule_off_list;

  p | prolong_type;
  p | restrict_type;

  p | stopping_cycle;
  p | stopping_time;
  p | stopping_interval;

  p | testing_cycle_final;
  p | testing_time_final;

  p | timestep_type;

}

//----------------------------------------------------------------------

void Config::read(Parameters * parameters) throw()
{
  TRACE("BEGIN Config::read()");

  //--------------------------------------------------
  // Boundary
  //--------------------------------------------------

  if (parameters->type("Boundary:list") != parameter_list) {

    num_boundary = 1;
    int ib = 0;

    // default name if not explicit

    boundary_list[ib] = "boundary";
    boundary_type[ib] = parameters->value_string("Boundary:type","unknown");

    std::string axis_str = parameters->value_string("Boundary:axis","all");
    if      (axis_str == "all") { boundary_axis[ib] = axis_all; }
    else if (axis_str == "x")   { boundary_axis[ib] = axis_x; }
    else if (axis_str == "y")   { boundary_axis[ib] = axis_y; }
    else if (axis_str == "z")   { boundary_axis[ib] = axis_z; }
    else {
      ERROR1 ("Config::read()", "Unknown Boundary:axis %s",
	      axis_str.c_str());
    }

    std::string face_str = parameters->value_string("Boundary:face","all");
    if      (face_str == "all")   { boundary_face[ib] = face_all; }
    else if (face_str == "lower") { boundary_face[ib] = face_lower; }
    else if (face_str == "upper") { boundary_face[ib] = face_upper; }
    else {
      ERROR1 ("Config::read()", "Unknown Boundary:face %s",
	      face_str.c_str());
    }

    boundary_mask[ib] = (parameters->type("Boundary:mask") 
		     == parameter_logical_expr);

    std::string param_str = "Boundary:field_list";
    int field_list_type = parameters->type(param_str);
    if (field_list_type == parameter_list) {
      const int n = parameters->list_length(param_str);
      boundary_field_list[ib].resize(n);
      for (int index=0; index<n; index++) {
	boundary_field_list[ib][index] = parameters->list_value_string 
	  (index,param_str);
      }
    } else if (field_list_type == parameter_string) {
      boundary_field_list[ib].resize(1);
      boundary_field_list[ib][0] = parameters->value_string(param_str);
    } else if (field_list_type != parameter_unknown) {
      ERROR2 ("Config::read()", "Incorrect parameter type %d for %s",
	      field_list_type,param_str.c_str());
    }
	

  } else {

    num_boundary = parameters->list_length("Boundary:list");

    ASSERT2 ("Config::read()", 
	     "Too many boundary conditions: %d must not be more than %d",
	     num_boundary,MAX_BOUNDARY,
	     num_boundary <= MAX_BOUNDARY);

    for (int ib=0; ib<num_boundary; ib++) {

      boundary_list[ib] = parameters->list_value_string
	(ib,"Boundary:list","unknown");

      std::string prefix = "Boundary:" + boundary_list[ib] + ":";
      boundary_type[ib] = parameters->value_string(prefix+"type","unknown");
      std::string axis_str = parameters->value_string(prefix+"axis","all");
      if      (axis_str == "all") { boundary_axis[ib] = axis_all; }
      else if (axis_str == "x")   { boundary_axis[ib] = axis_x; }
      else if (axis_str == "y")   { boundary_axis[ib] = axis_y; }
      else if (axis_str == "z")   { boundary_axis[ib] = axis_z; }
      else {
	ERROR2 ("Config::read()", "Unknown %s %s",
		(prefix+"axis").c_str(),axis_str.c_str());
      }

      std::string face_str = parameters->value_string(prefix+"face","all");
      if      (face_str == "all")   { boundary_face[ib] = face_all; }
      else if (face_str == "lower") { boundary_face[ib] = face_lower; }
      else if (face_str == "upper") { boundary_face[ib] = face_upper; }
      else {
	ERROR2 ("Config::read()", "Unknown %s %s",
		(prefix+"face").c_str(),face_str.c_str());
      }

      boundary_mask[ib] = (parameters->type(prefix+"mask") 
			  == parameter_logical_expr);

      std::string param_str = prefix+"field_list";
      int field_list_type = parameters->type(param_str);
      if (field_list_type == parameter_list) {
	const int n = parameters->list_length(param_str);
	boundary_field_list[ib].resize(n);
	for (int index=0; index<n; index++) {
	  boundary_field_list[ib][index] = parameters->list_value_string 
	    (index,param_str);
	}
      } else if (field_list_type == parameter_string) {
	boundary_field_list[ib].resize(1);
	boundary_field_list[ib][0] = parameters->value_string(param_str);
      } else if (field_list_type != parameter_unknown) {
	ERROR2 ("Config::read()", "Incorrect parameter type %d for %s",
		field_list_type,param_str.c_str());
      }
	
    }

  }

  //--------------------------------------------------
  // Domain
  //--------------------------------------------------

  for (int i=0; i<3; i++)  {
    domain_lower[i] = parameters->list_value_float(i, "Domain:lower", 0.0);
    domain_upper[i] = parameters->list_value_float(i, "Domain:upper", 0.0);
  }

  //--------------------------------------------------
  // Field
  //--------------------------------------------------

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
  // face_rank = rank - ( |ix| + |iy| + |iz| )
  //
  // face_rank               0     1     2         0      1
  //                    3D corner edge  face  2D corner edge
  // refresh_rank == 0:      x     x     x         x      x  >= 0D faces
  // refresh_rank == 1             x     x                x  >= 1D faces
  // refresh_rank == 2                   x                   >= 2D faces
  //
  // refresh if (face_rank >= rank)

  field_refresh_rank = parameters->value_integer ("Field:refresh:rank",0);

  prolong_type   = parameters->value_string ("Field:prolong","linear");

  restrict_type  = parameters->value_string ("Field:restrict","linear");

  //--------------------------------------------------
  // Initial
  //--------------------------------------------------

  TRACE("Parameters: Initial");
  initial_type  = parameters->value_string("Initial:type","default");
  initial_cycle = parameters->value_integer("Initial:cycle",0);
  initial_time  = parameters->value_float  ("Initial:time",0.0);
  initial_max_level = parameters->value_integer("Initial:max_level",0);

  //  initial_name;

  //  initial_value

  //--------------------------------------------------
  // Memory
  //--------------------------------------------------

  memory_active = parameters->value_logical("Memory:active",true);

  //--------------------------------------------------
  // Mesh
  //--------------------------------------------------

  TRACE("Parameters: Mesh");
  mesh_root_rank = parameters->value_integer("Mesh:root_rank",0);

  if (mesh_root_rank < 2) field_ghosts[1] = 0;
  if (mesh_root_rank < 3) field_ghosts[2] = 0;
  
  //--------------------------------------------------

  mesh_root_blocks[0] = parameters->list_value_integer(0,"Mesh:root_blocks",1);
  mesh_root_blocks[1] = parameters->list_value_integer(1,"Mesh:root_blocks",1);
  mesh_root_blocks[2] = parameters->list_value_integer(2,"Mesh:root_blocks",1);

  //--------------------------------------------------

  mesh_root_size[0] = parameters->list_value_integer(0,"Mesh:root_size",1);
  mesh_root_size[1] = parameters->list_value_integer(1,"Mesh:root_size",1);
  mesh_root_size[2] = parameters->list_value_integer(2,"Mesh:root_size",1);

  //--------------------------------------------------

  mesh_max_level = parameters->value_integer("Mesh:max_level",0);

  //--------------------------------------------------

  mesh_balance   = parameters->value_logical("Mesh:adapt:balance",true);

  //--------------------------------------------------

  int num_adapt_type = parameters->list_length("Mesh:adapt:type");
  mesh_adapt_type.resize(num_adapt_type);
  for (int i=0; i<num_adapt_type; i++) {
    mesh_adapt_type[i] = parameters->list_value_string(i,"Mesh:adapt:type");
  }

  //--------------------------------------------------

  control_sync_adapt_enter = parameters->value_string
    ("Control:sync:adapt_enter","array");
  control_sync_adapt_called = parameters->value_string
    ("Control:sync:adapt_called","neighbor");
  control_sync_adapt_end = parameters->value_string
    ("Control:sync:adapt_end","quiescence");
  control_sync_adapt_next = parameters->value_string
    ("Control:sync:adapt_next","quiescence");
  control_sync_adapt_exit = parameters->value_string
    ("Control:sync:adapt_exit","none");
  control_sync_compute_enter = parameters->value_string
    ("Control:sync:compute_enter","none");
  control_sync_compute_exit = parameters->value_string
    ("Control:sync:compute_exit","none");
  control_sync_exit = parameters->value_string
    ("Control:sync:exit = parameters","contribute");
  control_sync_output_enter = parameters->value_string
    ("Control:sync:output_enter","array");
  control_sync_output_exit = parameters->value_string
    ("Control:sync:output_exit","none");
  control_sync_refresh_enter = parameters->value_string
    ("Control:sync:refresh_enter","array");
  control_sync_refresh_exit = parameters->value_string
    ("Control:sync:refresh_exit","contribute");
  control_sync_stopping_enter = parameters->value_string
    ("Control:sync:stopping_enter","none");
  control_sync_stopping_exit = parameters->value_string
    ("Control:sync:stopping_exit","none");

  //--------------------------------------------------

  int num_adapt_fields = parameters->list_length("Mesh:adapt:fields");

  mesh_adapt_fields.resize(num_adapt_fields);

  for (int i=0; i<num_adapt_fields; i++) {
    mesh_adapt_fields[i] = parameters->list_value_string
      (i,"Mesh:adapt:fields");
  }

  //--------------------------------------------------

  mesh_adapt_slope_min_refine = 
    parameters->value_float ("Mesh:adapt:slope_min_refine",0.3);

  mesh_adapt_slope_max_coarsen = 
    parameters->value_float ("Mesh:adapt:slope_max_coarsen",0.15);

  //--------------------------------------------------

  // This parameter is typically computed ("internal" parameter in Enzo)
  mesh_adapt_mass_min = 
    parameters->value_float ("Mesh:adapt:mass_min",-1.0);

  //--------------------------------------------------

  mesh_adapt_mass_level_exponent = 
    parameters->value_float ("Mesh:adapt:mass_level_exponent",0.0);
  //--------------------------------------------------

  mesh_adapt_mass_min_overdensity = 
    parameters->value_float ("Mesh:adapt:mass_min_overdensity",1.5);

  mesh_adapt_balance = 
    parameters->value_logical ("Mesh:adapt:balance",true);

  mesh_adapt_interval = 
    parameters->value_integer ("Mesh:adapt:interval",1);

  //--------------------------------------------------
  // Method
  //--------------------------------------------------
 
  TRACE("Parameters: Method");

  int size = parameters->list_length("Method:sequence");
  method_sequence.resize(size);
  for (int i=0; i<size; i++) {
    method_sequence[i] = parameters->list_value_string(i,"Method:sequence");
  }

  //--------------------------------------------------
  // Monitor
  //--------------------------------------------------

  monitor_debug = parameters->value_logical("Monitor:debug",false);

  //--------------------------------------------------
  // Output
  //--------------------------------------------------

  parameters->group_set(0,"Output");

  num_file_groups = parameters->list_length("file_groups");

  ASSERT2 ("Config::read","Number of file groups %d exceeds MAX_FILE_GROUPS %d",
	   num_file_groups, MAX_FILE_GROUPS, num_file_groups <= MAX_FILE_GROUPS);

  parameters->group_set(0,"Output");


  output_file_groups.resize(num_file_groups);

  for (int index=0; index<num_file_groups; index++) {

    TRACE1 ("index = %d",index);

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
	TRACE3("output_dir[%d][%d] = %s",index,i,output_dir[index][i].c_str());
      }
    }

    TRACE1("index = %d",index);
    if (parameters->type("name") == parameter_string) {
      TRACE0;
      output_name[index].resize(1);
      output_name[index][0] = parameters->value_string("name","");
    } else if (parameters->type("name") == parameter_list) {
      int size = parameters->list_length("name");
      TRACE1("size = %d",size);
      if (size > 0) output_name[index].resize(size);
      for (int i=0; i<size; i++) {
	output_name[index][i] = parameters->list_value_string(i,"name","");
      }
    }

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


    read_schedule_(parameters, output_file_groups[index],
		   &output_schedule_type[index],
		   &output_schedule_var[index],
		   &output_schedule_start[index],
		   &output_schedule_stop[index],
		   &output_schedule_step[index],
		   output_schedule_list[index]
		   );

    // Image 

    TRACE2 ("output_type[%d] = %s",index,output_type[index].c_str());

    if (output_type[index] == "image") {

      WARNING1 ("Config::read()",
		"output_image_axis[%d] set to z",index);

      output_image_axis[index] = "z";

      if (parameters->type("axis") != parameter_unknown) {
	std::string axis = parameters->value_string("axis");
	ASSERT2("Problem::initialize_output",
		"Output %s axis %d must be \"x\", \"y\", or \"z\"",
		output_file_groups[index].c_str(), axis.c_str(),
		axis=="x" || axis=="y" || axis=="z");
      } 

      output_image_block_size[index] = 
	parameters->value_integer("image_block_size",1);

      output_image_type[index] = parameters->value_string("image_type","data");

      output_image_log[index] = parameters->value_logical("image_log",false);

      output_image_mesh_color[index] = 
	parameters->value_string("image_mesh_color","level");

      output_image_size[index].resize(2);
      output_image_size[index][0] = 
	parameters->list_value_integer(0,"image_size",0);
      output_image_size[index][1] = 
	parameters->list_value_integer(1,"image_size",0);

      output_image_reduce_type[index] = 
	parameters->value_string("image_reduce_type","sum");

      output_image_face_rank[index] = 
	parameters->value_integer("image_face_rank",3);

      output_image_ghost[index] = 
	parameters->value_logical("image_ghost",false);

      output_image_specify_bounds[index] =
	parameters->value_logical("image_specify_bounds",false);
      output_image_min[index] =
	parameters->value_float("image_min",0.0);
      output_image_max[index] =
	parameters->value_float("image_max",0.0);

      if (parameters->type("colormap") == parameter_list) {
	int size = parameters->list_length("colormap");
	output_image_colormap[index].resize(size);
	for (int i=0; i<size; i++) {
	  output_image_colormap[index][i] = 
	    parameters->list_value_float(i,"colormap",0.0);
	}
      }

    }
  }  

  //--------------------------------------------------
  // Performance
  //--------------------------------------------------

  if (parameters->type("Performance:papi:counters") == parameter_list) {
    int length = parameters->list_length("Performance:papi:counters");
    performance_papi_counters.resize(length);
    for (int i=0; i<length; i++) {
      performance_papi_counters[i] = parameters->list_value_string
	(i,"Performance:papi:counters","");
      TRACE2("performance_papi_counters[%d] = %s",
	     i,performance_papi_counters[i].c_str());
    }
  }

  performance_name     = parameters->value_string ("Performance:name","");
  performance_stride   = parameters->value_integer("Performance:stride",1);
  performance_warnings = parameters->value_logical("Performance:warnings",true);

  //
  parameters->group_set(0,"Performance");
  parameters->group_set(1,"projections");

  parameters->group_set(2,"on");
  read_schedule_(parameters,"on",
		 &projections_schedule_on_type,
		 &projections_schedule_on_var,
		 &projections_schedule_on_start,
		 &projections_schedule_on_stop,
		 &projections_schedule_on_step,
		 projections_schedule_on_list);

  
  parameters->group_set(2,"off");
  read_schedule_(parameters,"off",
		 &projections_schedule_off_type,
		 &projections_schedule_off_var,
		 &projections_schedule_off_start,
		 &projections_schedule_off_stop,
		 &projections_schedule_off_step,
		 projections_schedule_off_list);

  //--------------------------------------------------
  // Stopping
  //--------------------------------------------------

  stopping_cycle = parameters->value_integer
    ( "Stopping:cycle" , std::numeric_limits<int>::max() );
  stopping_time  = parameters->value_float
    ( "Stopping:time" , std::numeric_limits<double>::max() );
  stopping_interval = parameters->value_integer
    ( "Stopping:interval" , 1);

  //--------------------------------------------------
  // Testing
  //--------------------------------------------------

  testing_cycle_final = parameters->value_integer("Testing:cycle_final",0);
  testing_time_final  = parameters->value_float  ("Testing:time_final", 0.0);

  //--------------------------------------------------
  // Timestep
  //--------------------------------------------------

  timestep_type = parameters->value_string("Timestep:type","default");

  TRACE("END   Config::read()");

}

//======================================================================

void Config::read_schedule_(Parameters * parameters,
			    const std::string group,
			    std::string * type,
			    std::string * var,
			    double * start,
			    double * stop,
			    double * step,
			    std::vector<double> & list)
{
  // ASSUMES schedule PARAMETER IN CURRENT GROUP

  int length = parameters->list_length("schedule");

  if (length == 0) {
    (*type) = "none";
    return;
  }

  ASSERT1("Config::read","Incorrect 'schedule' for parameter group %s",
	  group.c_str(), (length >= 3));

  (*var)  = parameters->list_value_string(0,"schedule");
  (*type) = parameters->list_value_string(1,"schedule");

  bool var_is_int = true;

  if ((*var) == "cycle") {
    var_is_int = true;
  } else if ((*var) == "time") {
    var_is_int = false;
  } else {
    ERROR2 ("Config::read",
	    "Schedule variable %s is not recognized for parameter group %s",
	    (*var).c_str(),group.c_str());
  }

  const int    max_int    = std::numeric_limits<int>::max();
  const double max_double = std::numeric_limits<double>::max();

  if ((*type) == "interval") {
    if (length == 3) { // step
      if (var_is_int) {
	(*start) = 0;
	(*step)  = parameters->list_value_integer(2,"schedule");
	(*stop)  = max_int;
      } else {
	(*start) = 0.0;
	(*step)  = parameters->list_value_float(2,"schedule");
	(*stop)  = max_double;
      }
    } else if (length == 4) { // start, step
      if (var_is_int) {
	(*start) = parameters->list_value_integer(2,"schedule");
	(*step)  = parameters->list_value_integer(3,"schedule");
	(*stop)  = max_int;
      } else {
	(*start) = parameters->list_value_float(2,"schedule");
	(*step)  = parameters->list_value_float(3,"schedule");
	(*stop)  = max_double;
      }
    } else if (length == 5) { // start, step, stop
      if (var_is_int) {
	(*start) = parameters->list_value_integer(2,"schedule");
	(*step)  = parameters->list_value_integer(3,"schedule");
	(*stop)  = parameters->list_value_integer(4,"schedule");
      } else {
	(*start) = parameters->list_value_float(2,"schedule");
	(*step)  = parameters->list_value_float(3,"schedule");
	(*stop)  = parameters->list_value_float(4,"schedule");
      }
    }
  } else if ((*type) == "list") {
    list.resize(length-2);
    for (int i=2; i<length; i++) {
      if (var_is_int) {
	list[i-2] = parameters->list_value_integer(i,"schedule");
      } else {
	list[i-2] = parameters->list_value_float(i,"schedule");
      }
    }
  } else {
    ERROR2 ("Config::read",
	    "Schedule type %s is not recognized for parameter group %s",
	    (*type).c_str(),group.c_str());
  }

}

//======================================================================

