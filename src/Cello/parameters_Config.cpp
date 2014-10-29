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
  // REVIEW: rev 3505

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  // Balance

  p | balance_interval;

  // Boundary

  p | num_boundary;
  PUParray(p,boundary_list,MAX_BOUNDARY);
  PUParray(p,boundary_type,MAX_BOUNDARY);
  PUParray(p,boundary_axis,MAX_BOUNDARY);
  PUParray(p,boundary_face,MAX_BOUNDARY);
  PUParray(p,boundary_mask,MAX_BOUNDARY);
  PUParray(p,boundary_field_list,MAX_BOUNDARY);

  // Domain

  PUParray(p,domain_lower,3);
  PUParray(p,domain_upper,3);

  // Field

  p | num_fields;
  p | field_alignment;
  PUParray(p,field_centering,3);
  p | field_courant;
  p | field_list;
  PUParray(p,field_ghosts,3);
  p | field_padding;
  p | field_precision;
  p | field_refresh_rank;
  p | field_prolong_type;
  p | field_restrict_type;
  p | field_group_list;

  // Initial

  p | initial_cycle;
  p | initial_type;
  p | initial_time;

  // Memory

  p | memory_active;

  // Mesh

  PUParray(p,mesh_root_blocks,3);
  p | mesh_root_rank;
  PUParray(p,mesh_root_size,3);
  p | mesh_max_level;
  p | mesh_adapt_interval;
  p | num_mesh;
  PUParray(p,mesh_list,MAX_MESH_GROUPS);
  PUParray(p,mesh_type,MAX_MESH_GROUPS);
  PUParray(p,mesh_field_list,MAX_MESH_GROUPS);
  PUParray(p,mesh_min_refine,MAX_MESH_GROUPS);
  PUParray(p,mesh_max_coarsen,MAX_MESH_GROUPS);
  PUParray(p,mesh_min_refine2,MAX_MESH_GROUPS);
  PUParray(p,mesh_max_coarsen2,MAX_MESH_GROUPS);
  PUParray(p,mesh_level_exponent,MAX_MESH_GROUPS);
  PUParray(p,mesh_refine_output,MAX_MESH_GROUPS);

  // Method

  p | method_list;

  // Monitor

  p | monitor_debug;

  // Output

  p | num_output;
  p | output_list;
  PUParray (p,output_type,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_axis,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_block_size,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_colormap,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_type,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_log,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_mesh_color,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_size,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_reduce_type,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_ghost,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_face_rank,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_specify_bounds,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_min,MAX_OUTPUT_GROUPS);
  PUParray (p,output_image_max,MAX_OUTPUT_GROUPS);
  PUParray (p,output_schedule_index,MAX_OUTPUT_GROUPS);
  PUParray (p,output_field_list,MAX_OUTPUT_GROUPS);
  PUParray (p,output_stride,MAX_OUTPUT_GROUPS);
  PUParray (p,output_name,MAX_OUTPUT_GROUPS);
  PUParray (p,output_dir,MAX_OUTPUT_GROUPS);
  p | index_schedule_;
  PUParray (p,output_schedule_type,MAX_OUTPUT_GROUPS);
  PUParray (p,output_schedule_var,MAX_OUTPUT_GROUPS);
  PUParray (p,output_schedule_start,MAX_OUTPUT_GROUPS);
  PUParray (p,output_schedule_stop,MAX_OUTPUT_GROUPS);
  PUParray (p,output_schedule_step,MAX_OUTPUT_GROUPS);
  PUParray (p,output_schedule_list,MAX_OUTPUT_GROUPS);

  // Performance

  p | performance_papi_counters;
  p | performance_name;
  p | performance_stride;
  p | performance_warnings;

  // Restart

  p | restart_file;

  // Stopping

  p | stopping_cycle;
  p | stopping_time;
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
  read_method_(p);
  read_monitor_(p);
  read_output_(p);
  read_performance_(p);
  read_restart_(p);
  read_stopping_(p);
  read_testing_(p);

  TRACE("END   Config::read()");

}

//----------------------------------------------------------------------

void Config::read_balance_ (Parameters * p) throw()
{

  //--------------------------------------------------
  // Boundary
  //--------------------------------------------------

  balance_interval = p->value_integer("Balance:interval",0);

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


  ASSERT2 ("Config::read()", 
	   "Too many boundary conditions: %d must not be more than %d",
	   num_boundary,MAX_BOUNDARY,
	   num_boundary <= MAX_BOUNDARY);

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

  ASSERT2 ("Config::read","Number of fields %d exceeds MAX_FIELDS %d",
	   num_fields, MAX_FIELDS, num_fields <= MAX_FIELDS);

  field_list.resize(num_fields);
  for (int i=0; i<num_fields; i++) {
    field_list[i] = p->list_value_string(i, "Field:list");
  }

  if (p->type("Field:ghosts") == parameter_integer) {
    field_ghosts[0] = p->value_integer("Field:ghosts",0);
    field_ghosts[1] = p->value_integer("Field:ghosts",0);
    field_ghosts[2] = p->value_integer("Field:ghosts",0);
  } else if (p->type("Field:ghosts") == parameter_list) {
    field_ghosts[0] = p->list_value_integer(0,"Field:ghosts",0);
    field_ghosts[1] = p->list_value_integer(1,"Field:ghosts",0);
    field_ghosts[2] = p->list_value_integer(2,"Field:ghosts",0);
  } else {
    field_ghosts[0] = 0;
    field_ghosts[1] = 0;
    field_ghosts[2] = 0;
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
      int num_groups = p->list_length(param);
      for (int index_group=0; index_group<num_groups; index_group++) {
	std::string group = p->list_value_string(index_group,param);
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

  for (int index_group = 0; index_group < num_groups; index_group++) {

    std::string group = p->list_value_string(index_group, "Group:list") ;

    param = std::string("Group:") + group + ":field_list";

    if (p->type(param) == parameter_list) {
      // field_list is a list
      for (int i=0; i<num_fields; i++) {
	std::string field = p->list_value_string(i,param);
	for (int index_field = 0; index_field < num_fields; index_field++) {
	  if (field_list[index_field] == field) {
	    field_group_list[index_field].push_back(group);
	  }
	}
      }
    } else if (p->type(param) == parameter_string) {
      // field_list is a string
      std::string field = p->value_string(param);
      for (int index_field = 0; index_field < num_fields; index_field++) {
	if (field_list[index_field] == field) {
	  field_group_list[index_field].push_back(group);
	}
      }
    }
    
  }

  field_courant = p->value_float  ("Field:courant",0.6);

  field_padding = p->value_integer("Field:padding",0);

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
  // face_rank = rank - ( |ix| + |iy| + |iz| )
  //
  // face_rank               0     1     2         0      1
  //                    3D corner edge  face  2D corner edge
  // refresh_rank == 0:      x     x     x         x      x  >= 0D faces
  // refresh_rank == 1             x     x                x  >= 1D faces
  // refresh_rank == 2                   x                   >= 2D faces
  //
  // refresh if (face_rank >= rank)

  field_refresh_rank = p->value_integer ("Field:refresh:rank",0);

  field_prolong_type   = p->value_string ("Field:prolong","linear");

  field_restrict_type  = p->value_string ("Field:restrict","linear");
  
}

//----------------------------------------------------------------------

void Config::read_initial_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Initial
  //--------------------------------------------------

  TRACE("Parameters: Initial");
  initial_cycle = p->value_integer("Initial:cycle",0);
  initial_type  = p->value_string ("Initial:type","value");
  initial_time  = p->value_float  ("Initial:time",0.0);

  //  initial_name;

  //  initial_value

}

//----------------------------------------------------------------------

void Config::read_memory_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Memory
  //--------------------------------------------------

  memory_active = p->value_logical("Memory:active",true);

}

//----------------------------------------------------------------------

void Config::read_mesh_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Mesh
  //--------------------------------------------------

  TRACE("Parameters: Mesh");
  mesh_root_rank = p->value_integer("Mesh:root_rank",0);

  if (mesh_root_rank < 2) field_ghosts[1] = 0;
  if (mesh_root_rank < 3) field_ghosts[2] = 0;
  
  //--------------------------------------------------

  mesh_root_blocks[0] = p->list_value_integer(0,"Mesh:root_blocks",1);
  mesh_root_blocks[1] = p->list_value_integer(1,"Mesh:root_blocks",1);
  mesh_root_blocks[2] = p->list_value_integer(2,"Mesh:root_blocks",1);

  ASSERT4 ("Config::read_mesh_()",
	   "Number of root blocks %d x %d x %d must not be less than number of processes %d",
	   mesh_root_blocks[0],mesh_root_blocks[1],mesh_root_blocks[2],CkNumPes(),
	   mesh_root_blocks[0]*mesh_root_blocks[1]*mesh_root_blocks[2] >= CkNumPes());

  //--------------------------------------------------

  mesh_root_size[0] = p->list_value_integer(0,"Mesh:root_size",1);
  mesh_root_size[1] = p->list_value_integer(1,"Mesh:root_size",1);
  mesh_root_size[2] = p->list_value_integer(2,"Mesh:root_size",1);

}

//----------------------------------------------------------------------

void Config::read_adapt_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Adapt
  //--------------------------------------------------

  mesh_adapt_interval = p->value ("Adapt:interval",1);

  //--------------------------------------------------

  mesh_max_level = p->value_integer("Adapt:max_level",0);

  //--------------------------------------------------

  num_mesh = p->list_length("Adapt:list");

  for (int ia=0; ia<num_mesh; ia++) {

    mesh_list[ia] = p->list_value_string (ia,"Adapt:list","unknown");

    std::string prefix = "Adapt:" + mesh_list[ia] + ":";

    mesh_type[ia] = p->value_string(prefix+"type","unknown");

    std::string param_str = prefix + "field_list";

    int type = p->type(param_str);

    if (type == parameter_list) {
      const int n = p->list_length(param_str);
      mesh_field_list[ia].resize(n);
      for (int index=0; index<n; index++) {
	mesh_field_list[ia][index] = p->value(index,param_str,"none");
      }
    } else if (type == parameter_string) {
      mesh_field_list[ia].resize(1);
      mesh_field_list[ia][0] = p->value(param_str,"none");
    } else if (type != parameter_unknown) {
      ERROR2 ("Config::read()", "Incorrect parameter type %d for %s",
	      type,param_str.c_str());
    }

    //--------------------------------------------------

    if (p->type(prefix + "min_refine") == parameter_float) {

      mesh_min_refine[ia]  = p->value (prefix + "min_refine",0.3);
      mesh_max_coarsen[ia] = p->value (prefix + "max_coarsen",
				       0.5*mesh_min_refine[ia]);

    } else if (p->type(prefix + "min_refine") == parameter_list) {

      mesh_min_refine[ia]  = p->list_value_float (0,prefix + "min_refine",0.3);
      mesh_max_coarsen[ia] = p->list_value_float (0,prefix + "max_coarsen",
						  0.5*mesh_min_refine[ia]);

      mesh_min_refine2[ia]  = p->list_value_float (0,prefix + "min_refine",0.3);
      mesh_max_coarsen2[ia] = p->list_value_float (1,prefix + "max_coarsen",
						   0.5*mesh_min_refine2[ia]);
    }

    mesh_level_exponent[ia] = p->value (prefix + "level_exponent",0.0);

    mesh_refine_output[ia] = p->value_string (prefix + "output","");

  }
}

//----------------------------------------------------------------------

void Config::read_method_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Method
  //--------------------------------------------------
 
  TRACE("Parameters: Method");

  int size = p->list_length("Method:list");
  method_list.resize(size);
  for (int i=0; i<size; i++) {
    method_list[i] = p->list_value_string(i,"Method:list");
  }

}

//----------------------------------------------------------------------

void Config::read_monitor_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Monitor
  //--------------------------------------------------

  monitor_debug = p->value_logical("Monitor:debug",false);

}

//----------------------------------------------------------------------

void Config::read_output_ (Parameters * p) throw()
{
  //--------------------------------------------------
  // Output
  //--------------------------------------------------

  p->group_set(0,"Output");

  num_output = p->list_length("list");

  ASSERT2 ("Config::read","Number of file groups %d exceeds MAX_OUTPUT_GROUPS %d",
	   num_output, MAX_OUTPUT_GROUPS, num_output <= MAX_OUTPUT_GROUPS);

  p->group_set(0,"Output");


  output_list.resize(num_output);

  for (int index=0; index<num_output; index++) {

    TRACE1 ("index = %d",index);

    output_list[index] = 
      p->list_value_string (index,"Output:list","unknown");

    p->group_set(1,output_list[index]);

    output_type[index] = p->value_string("type","unknown");

    if (output_type[index] == "unknown") {
      ERROR1("Config::read",
	     "Output:%s:type parameter is undefined",
	     output_list[index].c_str());
    }

    output_stride[index] = p->value_integer("stride",0);

    if (p->type("dir") == parameter_string) {
      output_dir[index].resize(1);
      output_dir[index][0] = p->value_string("dir","");
    } else if (p->type("dir") == parameter_list) {
      int size = p->list_length("dir");
      if (size > 0) output_dir[index].resize(size);
      for (int i=0; i<size; i++) {
	output_dir[index][i] = p->list_value_string(i,"dir","");
	TRACE3("output_dir[%d][%d] = %s",index,i,output_dir[index][i].c_str());
      }
    }

    TRACE1("index = %d",index);
    if (p->type("name") == parameter_string) {
      TRACE0;
      output_name[index].resize(1);
      output_name[index][0] = p->value_string("name","");
    } else if (p->type("name") == parameter_list) {
      int size = p->list_length("name");
      TRACE1("size = %d",size);
      if (size > 0) output_name[index].resize(size);
      for (int i=0; i<size; i++) {
	output_name[index][i] = p->list_value_string(i,"name","");
      }
    }

    if (p->type("field_list") == parameter_list) {
      int length = p->list_length("field_list");
      output_field_list[index].resize(length);
      for (int i=0; i<length; i++) {
	output_field_list[index][i] = p->list_value_string(i,"field_list","");
      }
    }
      
    p->group_push("schedule");
    output_schedule_index[index] = 
      read_schedule_(p, output_list[index]);
    p->group_pop();

    // Image 

    TRACE2 ("output_type[%d] = %s",index,output_type[index].c_str());

    if (output_type[index] == "image") {

      WARNING1 ("Config::read()",
		"output_image_axis[%d] set to z",index);

      output_image_axis[index] = "z";

      if (p->type("axis") != parameter_unknown) {
	std::string axis = p->value_string("axis");
	ASSERT2("Problem::initialize_output",
		"Output %s axis %d must be \"x\", \"y\", or \"z\"",
		output_list[index].c_str(), axis.c_str(),
		axis=="x" || axis=="y" || axis=="z");
      } 

      output_image_block_size[index] = 
	p->value_integer("image_block_size",1);

      output_image_type[index] = p->value_string("image_type","data");

      output_image_log[index] = p->value_logical("image_log",false);

      output_image_mesh_color[index] = 
	p->value_string("image_mesh_color","level");

      output_image_size[index].resize(2);
      output_image_size[index][0] = 
	p->list_value_integer(0,"image_size",0);
      output_image_size[index][1] = 
	p->list_value_integer(1,"image_size",0);

      output_image_reduce_type[index] = 
	p->value_string("image_reduce_type","sum");

      output_image_face_rank[index] = 
	p->value_integer("image_face_rank",3);

      output_image_ghost[index] = 
	p->value_logical("image_ghost",false);

      output_image_specify_bounds[index] =
	p->value_logical("image_specify_bounds",false);
      output_image_min[index] =
	p->value_float("image_min",0.0);
      output_image_max[index] =
	p->value_float("image_max",0.0);

      if (p->type("colormap") == parameter_list) {
	int size = p->list_length("colormap");
	output_image_colormap[index].resize(size);
	for (int i=0; i<size; i++) {
	  output_image_colormap[index][i] = 
	    p->list_value_float(i,"colormap",0.0);
	}
      }

    }
  }  

}

//----------------------------------------------------------------------

void Config::read_performance_ (Parameters * p) throw()
{
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

  performance_name     = p->value_string ("Performance:name","");
  performance_stride   = p->value_integer("Performance:stride",1);
  performance_warnings = p->value_logical("Performance:warnings",true);

}

//----------------------------------------------------------------------

void Config::read_restart_ (Parameters * p) throw()
{
  restart_file = p->value_string("Restart:file","");
}

//----------------------------------------------------------------------

void Config::read_stopping_ (Parameters * p) throw()
{
  stopping_cycle = p->value_integer
    ( "Stopping:cycle" , std::numeric_limits<int>::max() );
  stopping_time  = p->value_float
    ( "Stopping:time" , std::numeric_limits<double>::max() );
  stopping_interval = p->value_integer
    ( "Stopping:interval" , 1);
}

//----------------------------------------------------------------------

void Config::read_testing_ (Parameters * p) throw()
{
  testing_cycle_final = p->value_integer("Testing:cycle_final",0);
  testing_time_final  = p->value_float  ("Testing:time_final", 0.0);
  testing_time_tolerance = p->value_float  ("Testing:time_tolerance", 1e-6);
}

//======================================================================

int Config::read_schedule_(Parameters * p, const std::string group)
{
  int index = index_schedule_;
  ASSERT ("Config::read_schedule_()",
	  "number of schedule''s is greater than MAX_SCHEDULE",
	  index < MAX_SCHEDULE);

  std::string var = p->value_string("var","none");

  output_schedule_var[index] = var;

  bool var_is_int = true;

  // Get variable associated with the schedule 
  if      (output_schedule_var[index] == "cycle")    var_is_int = true;
  else if (output_schedule_var[index] == "time")     var_is_int = false;
  else if (output_schedule_var[index] == "seconds")  var_is_int = false;
  else {
    ERROR2 ("Config::read",
	    "Schedule variable %s is not recognized for parameter group %s",
	    output_schedule_var[index].c_str(),group.c_str());
  }

  // Determine the schedule type (interval or list)

  const bool type_is_interval = 
    ( (p->type("start") != parameter_unknown) ||
      (p->type("step") != parameter_unknown) ||
      (p->type("stop") != parameter_unknown));
  
  const bool type_is_list = (p->type("value") != parameter_unknown);

  if (type_is_interval && type_is_list) {
      ERROR1 ("Config::read",
	      "Schedule %s seems to be both an interval and a list",
	      (group).c_str());
  }

  if (!type_is_interval && !type_is_list) {
      ERROR1 ("Config::read",
	      "Schedule %s seems to be neither an interval nor a list",
	      (group).c_str());
  }

  if (type_is_interval) output_schedule_type[index] = "interval";
  if (type_is_list)     output_schedule_type[index] = "list";

  const int    max_int    = std::numeric_limits<int>::max();
  const double max_double = std::numeric_limits<double>::max();

  if (type_is_interval) {
    if (var_is_int) {
      output_schedule_start[index] = p->value("start",0);
      output_schedule_step[index]  = p->value("step",1);
      output_schedule_stop[index]  = p->value("stop",max_int);
    } else {
      output_schedule_start[index] = p->value("start",0.0);
      output_schedule_step[index]  = p->value("step",1.0);
      output_schedule_stop[index]  = p->value("stop",max_double);
    }
  } else if (type_is_list) {
    int n = p->list_length("value");
    if (n == 0) {
      ERROR1 ("Config::read",
	      "Schedule variable %s has length 0",
	      (group + ":list").c_str());
    }
    output_schedule_list[index].resize(n);
    if (var_is_int) {
      for (int i=0; i<n; i++) {
	output_schedule_list[index][i] = p->value(i,"value",0);
      }
    } else {
      for (int i=0; i<n; i++) {
	output_schedule_list[index][i] = p->value(i,"value",0.0);
      }
    }
  } else {
    ERROR2 ("Config::read",
	    "Schedule type %s is not recognized for parameter group %s",
	    output_schedule_type[index].c_str(),group.c_str());
  }

  return index_schedule_++;
}

//======================================================================

