// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Problem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-03
/// @brief    Implementation of the Problem container class

#include "problem.hpp"

//----------------------------------------------------------------------

Problem::Problem() throw()
  : boundary_(0),
    initial_(0),
    stopping_(0),
    timestep_(0)
{
}

//----------------------------------------------------------------------

Problem::~Problem() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Problem::initialize_boundary(Parameters * parameters) throw()
{
  //--------------------------------------------------
  // parameter: Boundary : type
  //--------------------------------------------------

  TRACE0;
  std::string type = parameters->value_string("Boundary:type","");

  boundary_ = create_boundary_(type,parameters);

  ASSERT1("Problem::initialize_boundary",
	  "Boundary type %s not recognized",
	  type.c_str(),
	  boundary_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_initial(Parameters * parameters,
				 GroupProcess * group_process) throw()
{
  TRACE0;
  //--------------------------------------------------
  // parameter: Initial : type
  //--------------------------------------------------

  std::string type = parameters->value_string("Initial:type","default");

  initial_ = create_initial_(type,parameters,group_process);

  ASSERT1("Problem::initialize_initial",
	  "Initial type %s not recognized",
	  type.c_str(),
	  initial_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_stopping(Parameters * parameters) throw()
{
  TRACE0;
  stopping_ = create_stopping_("default",parameters);

  ASSERT("Problem::initialize_stopping",
	  "Stopping object not successfully created",
	  stopping_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_timestep(Parameters * parameters) throw()
{
  TRACE0;
  //--------------------------------------------------
  // parameter: Timestep : type
  //--------------------------------------------------

  std::string type = parameters->value_string("Timestep:type","default");

  timestep_ = create_timestep_(type,parameters);

  ASSERT1("Problem::initialize_timestep",
	  "Timestep type %s not recognized",
	  type.c_str(),
	  timestep_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_output
(Parameters * parameters,
 FieldDescr * field_descr,
 GroupProcess * group_process,
 Hierarchy    * hierarchy,
 const Factory * factory) throw()
{
  TRACE0;
  // Create and initialize an Output object for each Output group

  //--------------------------------------------------
  // parameter: Output : file_groups
  //--------------------------------------------------

  //--------------------------------------------------
  parameters->group_set(0,"Output");
  //--------------------------------------------------

  int num_file_groups = parameters->list_length("file_groups");


  for (int index_file_group=0; index_file_group < num_file_groups; index_file_group++) {

    //--------------------------------------------------
    parameters->group_set(0,"Output");
    //--------------------------------------------------

    std::string file_group = parameters->list_value_string
      (index_file_group,"file_groups","unknown");

    //--------------------------------------------------
    parameters->group_set(1,file_group);

    //--------------------------------------------------
    // File type parameter
    //--------------------------------------------------

    //--------------------------------------------------
    // parameter: Output : <file_group> : type
    //--------------------------------------------------

    std::string type = parameters->value_string("type","unknown");

    // Error if Output::type is not defined
    if (type == "unknown") {
      ERROR1("Problem::initialize_output",
	     "Output:%s:type parameter is undefined",
	     file_group.c_str());
    }

    // Create output object

    Output * output = create_output_
      (type,parameters,group_process,hierarchy,factory);

    // Error if output type was not recognized
    if (output == NULL) {
      ERROR2("Problem::initialize_output",
	     "Unknown parameter type Output:%s:type = %s",
	     file_group.c_str(),
	     type.c_str());
    }

    //--------------------------------------------------
    // File name parameter
    //--------------------------------------------------

    std::string file_name = "";
    std::vector<std::string> file_args;

    //--------------------------------------------------
    // parameter: Output : <file_group> : name
    //--------------------------------------------------

    if (parameters->type("name") == parameter_string) {

      // Case 1: string e.g. name = "filename";

      file_name = parameters->value_string("name","");

     
    } else if (parameters->type("name") == parameter_list) {
      // Case 2: list e.g. name = ["filename-cycle-%0d.time-%5.3f", "cycle","time"]

      int list_length = parameters->list_length("name");

      // get file name
      if (list_length > 0) {
	file_name = parameters->list_value_string(0,"name","");
      }

      // get file args ("cycle", "time", etc.) to schedule
      for (int index = 1; index<list_length; index++) {
	file_args.push_back(parameters->list_value_string(index,"name",""));
      }

    } else {

      ERROR2("Problem::initialize_output",
	     "Bad type %d for 'Output : %s : name' parameter",
	     parameters->type("name"),file_group.c_str());

    }

    // Error check name specified

    ASSERT1("Problem::initialize_output",
	   "Output 'name' must be specified for file group %s",
	   file_group.c_str(),
	   file_name != "");

    // Set the output object file name and arguments

    output->set_filename (file_name,file_args);

    //--------------------------------------------------
    // Field_list parameter
    //--------------------------------------------------

    //--------------------------------------------------
    // parameter: Output : <file_group> : field_list
    //--------------------------------------------------

    if (parameters->type("field_list") != parameter_unknown) {

      // Set field list to specified field list

      if (parameters->type("field_list") == parameter_list) {

	ItFieldList * it_field = new ItFieldList;

	int length = parameters->list_length("field_list");
	for (int i=0; i<length; i++) {
	  std::string field_name = 
	    parameters->list_value_string(i,"field_list","");
	  int field_index = field_descr->field_id(field_name);
	  it_field->append(field_index);
	}
	
	output->set_it_field(it_field);

      } else {

	ERROR("Problem::initialize_output",
	      "Bad type for Output 'field_list' parameter");

      }

    } else {

      // field_list not defined: default to all fields
      int field_count = field_descr->field_count();
      ItFieldRange * it_field = new ItFieldRange(field_count);

      output->set_it_field(it_field);

    }

    //--------------------------------------------------
    // Scheduling parameters
    //--------------------------------------------------

    //--------------------------------------------------
    // parameter: Output : <file_group> : schedule
    //--------------------------------------------------

    // error check schedule parameter exists

    ASSERT("Problem::initialize_output",
	   "The 'schedule' parameter must be specified for all Output file groups",
	   parameters->type("schedule") != parameter_unknown);

    // get schedule variable ("cycle" or "time")

    bool var_cycle,var_time;
    var_cycle = (strcmp(parameters->list_value_string(0,"schedule"),"cycle")==0);
    var_time  = (strcmp(parameters->list_value_string(0,"schedule"),"time")==0);

    // error check variable name

    ASSERT("Problem::initialize_output",
	   "The first 'schedule' parameter list element must be 'cycle' or 'time'",
	   var_cycle || var_time);

    // get schedule type (interval or list)

    bool type_interval,type_list;

    type_interval = (strcmp(parameters->list_value_string(1,"schedule"),"interval")==0);
    type_list     = (strcmp(parameters->list_value_string(1,"schedule"),"list")==0);

    // error check schedule type

    ASSERT("Problem::initialize_output",
	   "The second 'schedule' parameter list element "
	   "must be 'interval' or 'list'",
	   type_interval || type_list);

    int len = parameters->list_length("schedule");

    if (var_cycle && type_interval) {

      // get cycle interval schedule

      const int max_int = std::numeric_limits<int>::max();

      // Set time limits

      int start = 0, stop = 0, step = 0;

      if (len == 3) { 

	// schedule = [ "cycle", "interval",  <step> ];

	start = 0;
	stop  = max_int;
	step  = parameters->list_value_integer(2,"schedule");

      } else if (len == 5) {

	// schedule = [ "cycle", "interval",  <start>, <stop>, <step> ]

	start = parameters->list_value_integer(2,"schedule");
	stop  = parameters->list_value_integer(3,"schedule");
	step  = parameters->list_value_integer(4,"schedule");

      } else {

	// error check schedule parameter list length

	ERROR("Problem::initialize_output",
	      "Output 'schedule' list parameter has wrong number of elements");

      }

      // Set cycle interval output schedule 

      output->schedule()->set_cycle_interval(start,step,stop);

    } else if (var_cycle && type_list) {

      // get cycle list schedule

      std::vector<int> list;

      for (int index = 2; index < len; index++) {

	int value = parameters->list_value_integer(index,"schedule");

	list.push_back (value);

	if (list.size() > 1) {

	  // error check monotonicity

	  ASSERT("Problem::initialize_output",
		 "Output 'schedule' parameter list values must be monotonically increasing",
		 (list[list.size()-2] < list[list.size()-1]));

	}
      }

      // Set cycle list output schedule 
      
      output->schedule()->set_cycle_list(list);

    } else if (var_time && type_interval) {

      // get time interval schedule

      const double max_double = std::numeric_limits<double>::max();

      // set time limits

      double start = 0, stop = -1, step = 1;

      if (len == 3) {

	// schedule = [ "time", "interval",  <step> ];

	start = 0;
	stop  = max_double;
	step  = parameters->list_value_float(2,"schedule");

      } else if (len == 5) {

	// schedule = [ "time", "interval",  <start>, <stop>, <step> ]

	start = parameters->list_value_float(2,"schedule");
	stop  = parameters->list_value_float(3,"schedule");
	step  = parameters->list_value_float(4,"schedule");

      } else {

	// error check schedule parameter list length

	ERROR("Problem::initialize_output",
	      "Output 'schedule' list parameter has wrong number of elements");

      }

      // Set time interval output schedule 

      output->schedule()->set_time_interval(start,step,stop);

    } else if (var_time && type_list) {

      // get time list schedule

      std::vector<double> list;

      for (int index = 2; index < len; index++) {

	double value = parameters->list_value_float(index,"schedule");

	list.push_back (value);

	if (list.size() > 1) {

	  // error check monotonicity

	  ASSERT("Problem::initialize_output",
		 "Output 'schedule' parameter list values must be "
		 "monotonically increasing",
		 (list[list.size()-2] < list[list.size()-1]));
	}
      }

      // Set time list output schedule 

      output->schedule()->set_time_list(list);

    }

    //--------------------------------------------------
    // Image parameters
    //--------------------------------------------------

    OutputImage * output_image = dynamic_cast<OutputImage *> (output);

    // WARNING("Problem::initialize_output()",
    // 	    "Temporarily setting output_image ghosts");
    //    output_image->set_ghosts(3,3,0);
    if (output != NULL) {

      // Parse image-specific parameters

      parameter_enum type = parameter_unknown;

      // position parameter

      type = parameters->type("position");

      if ( type != parameter_unknown) {

	// error check position type

	ASSERT1("Problem::initialize_output",
	       "Output %s position must be a float",
		file_group.c_str(),
		type == parameter_float);

	INCOMPLETE("Problem::initialize_output_(): image position not implemented");

      }

      // axis parameter

      type = parameters->type("axis");

      //--------------------------------------------------
      // parameter: Output : <file_group> : axis
      //--------------------------------------------------

      if (type != parameter_unknown) {

	// error check position type

	ASSERT1("Problem::initialize_output",
		"Output %s axis must be a string",
		file_group.c_str(), type == parameter_string);

	std::string axis = parameters->value_string("axis");

	// set the output axis

	if (axis == "x") output_image->set_axis(axis_x);
	if (axis == "y") output_image->set_axis(axis_y);
	if (axis == "z") output_image->set_axis(axis_z);

	// error check axis

	ASSERT2("Problem::initialize_output",
		"Output %s axis %d must be \"x\", \"y\", or \"z\"",
		file_group.c_str(), axis.c_str(),
		axis=="x" || axis=="y" || axis=="z");


      }

      //--------------------------------------------------
      // parameter: Output : <file_group> : colormap
      //--------------------------------------------------

      type = parameters->type("colormap");

      if (type != parameter_unknown) {

	// error check colormap list type

	ASSERT1("Problem::initialize_output",
		"Output %s colormap must be a list",
		file_group.c_str(), type == parameter_list);

	int n = parameters->list_length("colormap");

	// error check colormap list length

	ASSERT1("Problem::initialize_output",
		"Output %s colormap list length must be divisible by 3",
		file_group.c_str(), n % 3 == 0);

	n /= 3;

	// allocate arrays

	double * r = new double [n];
	double * g = new double [n];
	double * b = new double [n];


	for (int i=0; i<n; i++) {


	  int ir=3*i+0;
	  int ig=3*i+1;
	  int ib=3*i+2;

	  // error check colormap value types

	  ASSERT1("Problem::initialize_output",
		  "Output %s colormap list must only contain floats",
		  file_group.c_str(), 
		  ((parameters->list_type(ir,"colormap") == parameter_float) &&
		   (parameters->list_type(ig,"colormap") == parameter_float) &&
		   (parameters->list_type(ib,"colormap") == parameter_float)));

	  // get next colormap r[i] g[i], b[i]
	  r[i] = parameters->list_value_float (ir, "colormap",0.0);
	  g[i] = parameters->list_value_float (ig, "colormap",0.0);
	  b[i] = parameters->list_value_float (ib, "colormap",0.0);

	}

	// set the colormap

	output_image->set_colormap(n,r,g,b);

	// deallocate arrays

	delete r;
	delete g;
	delete b;

      }

      //--------------------------------------------------
      // parameter: Output : <file_group> : colormap_alpha
      //--------------------------------------------------

      type = parameters->type("colormap_alpha");

      if (type != parameter_unknown) {

	// error check colormap_alpha list type

	ASSERT1("Problem::initialize_output",
		"Output %s colormap_alpha must be a list",
		file_group.c_str(), type == parameter_list);

	int n = parameters->list_length("colormap_alpha");

	// error check colormap_alpha list length

	ASSERT1("Problem::initialize_output",
		"Output %s colormap_alpha list length must be divisible by 4",
		file_group.c_str(), n % 4 == 0);

	n /= 4;

	// allocate arrays

	double * r = new double [n];
	double * g = new double [n];
	double * b = new double [n];
	double * a = new double [n];

	for (int i=0; i<n; i++) {


	  int ir=4*i+0;
	  int ig=4*i+1;
	  int ib=4*i+2;
	  int ia=4*i+3;

	  // error check colormap_alpha value types

	  ASSERT1("Problem::initialize_output",
		  "Output %s colormap_alpha list must only contain floats",
		  file_group.c_str(), 
		  ((parameters->list_type(ir,"colormap_alpha") == parameter_float) &&
		   (parameters->list_type(ig,"colormap_alpha") == parameter_float) &&
		   (parameters->list_type(ib,"colormap_alpha") == parameter_float) &&
		   (parameters->list_type(ia,"colormap_alpha") == parameter_float)));

	  // get next colormap r[i] g[i], b[i]
	  r[i] = parameters->list_value_float (ir, "colormap_alpha",0.0);
	  g[i] = parameters->list_value_float (ig, "colormap_alpha",0.0);
	  b[i] = parameters->list_value_float (ib, "colormap_alpha",0.0);
	  a[i] = parameters->list_value_float (ia, "colormap_alpha",0.0);

	}

	// set the colormap

	output_image->set_colormap(n,r,g,b,a);

	// deallocate arrays

	delete r;
	delete g;
	delete b;
	delete a;

      }


    }

    // Initialize index of Output object in Simulation for CHARM
    // Output read()/write() callback by Hierarchy, Patch or Block

#ifdef CONFIG_USE_CHARM
    output->set_index_charm(output_list_.size());
#endif

    // Add the initialized Output object to the Simulation's list of
    // output objects

    output_list_.push_back(output); 


  } // (for index_file_group)

}

//----------------------------------------------------------------------

void Problem::initialize_method(Parameters * parameters) throw()
{
  TRACE0;

  //--------------------------------------------------
  parameters->group_set(0,"Method");
  //--------------------------------------------------

  int count = parameters->list_length("sequence");

  for (int i=0; i<count; i++) {

    //--------------------------------------------------
    // parameter: Method : sequence
    //--------------------------------------------------

    std::string name = parameters->list_value_string(i,"sequence");

    Method * method = create_method_(name,parameters);

    if (method) {

      method_list_.push_back(method); 

    } else {
      ERROR1("Problem::initialize_method",
	     "Unknown Method %s",name.c_str());
    }
  }
}

//======================================================================

void Problem::deallocate_() throw()
{
  delete boundary_;      boundary_ = 0;
  delete initial_;       initial_ = 0;
  delete stopping_;      stopping_ = 0;
  delete timestep_;      timestep_ = 0;
  for (size_t i=0; i<output_list_.size(); i++) {
    delete output_list_[i];    output_list_[i] = 0;
  }
  for (size_t i=0; i<method_list_.size(); i++) {
    delete method_list_[i];    method_list_[i] = 0;
  }
}

//----------------------------------------------------------------------

Boundary * Problem::create_boundary_
(
 std::string  name,
 Parameters * parameters
 ) throw ()
{
  // No default Boundary object
  return NULL;
}

//----------------------------------------------------------------------

Initial * Problem::create_initial_
(
 std::string  type,
 Parameters * parameters,
 GroupProcess * group_process
 ) throw ()
{ 
  //--------------------------------------------------
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

  int    init_cycle  = parameters->value_integer ("Initial:cycle",0);
  double init_time   = parameters->value_float   ("Initial:time",0.0);


  if (type == "file" || type == "restart") {
    return new InitialFile(parameters,group_process,init_cycle,init_time);;
  } else if (type == "default") {
    return new InitialDefault(parameters,init_cycle,init_time);
  }
  return NULL;
}

//----------------------------------------------------------------------

Stopping * Problem::create_stopping_ 
(
 std::string  type,
 Parameters * parameters
 ) throw ()
/// @param type   Type of the stopping method to create (ignored)
/// @param stop_cycle  Stopping cycle
/// @param stop_time  Stopping time
{
  //--------------------------------------------------
  // parameter: Stopping : cycle
  // parameter: Stopping : time
  //--------------------------------------------------

  int    stop_cycle = parameters->value_integer
    ( "Stopping:cycle" , std::numeric_limits<int>::max() );
  double stop_time  = parameters->value_float
    ( "Stopping:time" , std::numeric_limits<double>::max() );

  // Return default stopping criteria object

  return new Stopping(stop_cycle,stop_time);
}

//----------------------------------------------------------------------

Timestep * Problem::create_timestep_ 
(
 std::string  type,
 Parameters * parameters
 ) throw ()
{ 
  // No default timestep
  return NULL;
}

//----------------------------------------------------------------------

Method * Problem::create_method_ 
(
 std::string  type,
 Parameters * parameters
 ) throw ()
{
  // No default method
  return NULL;
}


//----------------------------------------------------------------------

Output * Problem::create_output_ 
(
 std::string    type,
 Parameters *   parameters,
 GroupProcess * group_process,
 Hierarchy    * hierarchy,
 const Factory * factory
 ) throw ()
/// @param type          Type of Output object to create
/// @param group_process Image output needs group process size
/// @param hierarchy     Image output needs image size (currently patch(0) size)
{ 
  Output * output = NULL;
  if (type == "image") {
    int np = group_process->size();
    int nx,ny,nz;
    hierarchy->patch(0)->size (&nx, &ny, &nz);
    // NOTE: assumes cube for non-z axis images
    output = new OutputImage (factory,np,nx,ny);
  } else if (type == "data") {
    output = new OutputData (factory);
  } else if (type == "restart") {
    output = new OutputRestart (factory,parameters);
  }
  return output;
}

