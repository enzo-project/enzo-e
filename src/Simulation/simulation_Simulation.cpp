// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "main.hpp"

#include "simulation.hpp"
#include "simulation_charm.hpp"

Simulation::Simulation
(
 const char *   parameter_file,
#ifdef CONFIG_USE_CHARM
 int n,
 CProxy_BlockReduce proxy_block_reduce,
#else
 GroupProcess * group_process,
#endif
 int            index
 ) throw()
/// Initialize the Simulation object
  : factory_(0),
    parameters_(0),
    parameter_file_(parameter_file),
#ifdef CONFIG_USE_CHARM
    group_process_(GroupProcess::create()),
#else
    group_process_(group_process),
#endif
    dimension_(0),
    cycle_(0),
    time_(0.0),
    dt_(0),
    stop_(false),
    temp_update_all_(true),
    temp_update_full_(true),
    index_(index),
    performance_(0),
    monitor_(0),
    hierarchy_(0),
    field_descr_(0),
    stopping_(0),
    timestep_(0),
    initial_(0),
    boundary_(0),
#ifdef CONFIG_USE_CHARM
    proxy_block_reduce_(proxy_block_reduce),
    index_output_(0),
#endif
    output_list_(),
    method_list_()
{
  performance_ = new Performance;
#ifdef CONFIG_USE_CHARM
  monitor_ = new Monitor;
#else
  monitor_ = Monitor::instance();
#endif
  parameters_ = new Parameters(parameter_file,monitor_);
}

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Simulation::initialize() throw()
{
  performance_->start();

  // Initialize parameters

  initialize_simulation_();

  // Initialize simulation components

  initialize_data_();
  initialize_hierarchy_();
  initialize_stopping_();
  initialize_timestep_();
  initialize_initial_();
  initialize_boundary_();
  initialize_output_();
  initialize_method_();
  initialize_parallel_();

  char parameter_file[40];
  snprintf (parameter_file,40,"%s.out",parameter_file_.c_str());
  parameters_->write(parameter_file);
}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  performance_->stop();
  deallocate_();
}

//======================================================================

void Simulation::initialize_simulation_() throw()
{

  //--------------------------------------------------
  // parameter: Physics::dimensions
  //--------------------------------------------------

  dimension_ = parameters_->value_integer("Physics:dimensions",0);

}

//----------------------------------------------------------------------

void Simulation::initialize_data_() throw()
{

  /* MEMORY LEAK */
  field_descr_ = new FieldDescr;

  //--------------------------------------------------
  // parameter: Field::fields
  //--------------------------------------------------

  // Add data fields

  int i;
  for (i=0; i<parameters_->list_length("Field:fields"); i++) {
    field_descr_->insert_field
      (parameters_->list_value_string(i, "Field:fields"));
  }

  // Define default ghost zone depth for all fields, default value of 1

  //--------------------------------------------------
  // parameter: Field::ghosts
  //--------------------------------------------------

  int gx = 1;
  int gy = 1;
  int gz = 1;

  if (parameters_->type("Field:ghosts") == parameter_integer) {
    gx = gy = gz = parameters_->value_integer("Field:ghosts",1);
  } else if (parameters_->type("Field:ghosts") == parameter_list) {
    gx = parameters_->list_value_integer(0,"Field:ghosts",1);
    gy = parameters_->list_value_integer(1,"Field:ghosts",1);
    gz = parameters_->list_value_integer(2,"Field:ghosts",1);
  }

  if (dimension_ < 2) gy = 0;
  if (dimension_ < 3) gz = 0;

  for (i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_ghosts(i,gx,gy,gz);
    int gx2,gy2,gz2;
    field_descr_->ghosts(i,&gx2,&gy2,&gz2);
  }

  // Set precision

  //--------------------------------------------------
  // parameter: Field::precision
  //--------------------------------------------------

  std::string precision_str = 
    parameters_->value_string("Field:precision","default");

  precision_enum precision = precision_default;

  if (precision_str == "single")
    precision = precision_single;
  else if (precision_str == "double")
    precision = precision_double;
  else if (precision_str == "quadruple")
    precision = precision_quadruple;

  for (i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_precision(i,precision);
  }
}

//----------------------------------------------------------------------

void Simulation::initialize_hierarchy_() throw()
{

  ASSERT("Simulation::initialize_hierarchy_",
	 "data must be initialized before hierarchy",
	 field_descr_ != NULL);

  //----------------------------------------------------------------------
  // Create and initialize Hierarchy
  //----------------------------------------------------------------------

  hierarchy_ = factory().create_hierarchy();

  // Domain extents

  //--------------------------------------------------
  // parameter: Domain::lower
  // parameter: Domain::upper
  //--------------------------------------------------

  ASSERT ("Simulation::initialize_simulation_",
	  "Parameter Domain:lower list length must match Physics::dimension",
	  (parameters_->list_length("Domain:lower") == dimension_));

  ASSERT ("Simulation::initialize_simulation_",
	  "Parameter Domain:upper list length must match Physics::dimension",
	  (parameters_->list_length("Domain:upper") ==  dimension_));

  double lower[3];
  double upper[3];

  for (int i=0; i<3; i++) {
    lower[i] = parameters_->list_value_float(i, "Domain:lower", 0.0);
    upper[i] = parameters_->list_value_float(i, "Domain:upper", 0.0);
    ASSERT ("Simulation::initialize_simulation_",
	    "Domain:lower may not be greater than Domain:upper",
	    lower[i] <= upper[i]);
  }

  hierarchy_->set_lower(lower[0], lower[1], lower[2]);
  hierarchy_->set_upper(upper[0], upper[1], upper[2]);

  //----------------------------------------------------------------------
  // Create and initialize root Patch in Hierarchy
  //----------------------------------------------------------------------

  //--------------------------------------------------
  // parameter: Mesh::root_size
  // parameter: Mesh::root_blocks
  //--------------------------------------------------

  int root_size[3];

  root_size[0] = parameters_->list_value_integer(0,"Mesh:root_size",1);
  root_size[1] = parameters_->list_value_integer(1,"Mesh:root_size",1);
  root_size[2] = parameters_->list_value_integer(2,"Mesh:root_size",1);

  int root_blocks[3];

  root_blocks[0] = parameters_->list_value_integer(0,"Mesh:root_blocks",1);
  root_blocks[1] = parameters_->list_value_integer(1,"Mesh:root_blocks",1);
  root_blocks[2] = parameters_->list_value_integer(2,"Mesh:root_blocks",1);

  hierarchy_->create_root_patch
    (group_process_,
     field_descr_,
     root_size[0],root_size[1],root_size[2],
     root_blocks[0],root_blocks[1],root_blocks[2]);

}

//----------------------------------------------------------------------

void Simulation::initialize_stopping_() throw()
{
  stopping_ = create_stopping_("ignored");
}

//----------------------------------------------------------------------

void Simulation::initialize_timestep_() throw()
{
  timestep_ = create_timestep_("ignored");
}

//----------------------------------------------------------------------

void Simulation::initialize_initial_() throw()
{
  //--------------------------------------------------
  // parameter: Initial::code
  //--------------------------------------------------

  std::string name = parameters_->value_string("Initial:code","default");

  initial_ = create_initial_(name);

  if (initial_ == 0) {
    initial_ = new InitialDefault(parameters_);
  }
}

//----------------------------------------------------------------------

void Simulation::initialize_boundary_() throw()
{
  //--------------------------------------------------
  // parameter: Boundary::name
  //--------------------------------------------------

  std::string name = parameters_->value_string("Boundary:name","");
  boundary_ = create_boundary_(name);
}

//----------------------------------------------------------------------

void Simulation::initialize_output_() throw()
{
  // Create and initialize an Output object for each Output group

  //--------------------------------------------------
  parameters_->group_set(0,"Output");
  //--------------------------------------------------

  int num_file_groups = parameters_->list_length("file_groups");

  for (int index_file_group=0; index_file_group < num_file_groups; index_file_group++) {

    //--------------------------------------------------
    parameters_->group_set(0,"Output");
    //--------------------------------------------------

    std::string file_group = parameters_->list_value_string
      (index_file_group,"file_groups","unknown");

    //--------------------------------------------------
    parameters_->group_set(1,file_group);
    //--------------------------------------------------
    // parameter: Output:<file_group>:type
    // parameter: Output:<file_group>:file_name
    // parameter: Output:<file_group>:field_list
    // parameter: Output:<file_group>:schedule
    //--------------------------------------------------

    //--------------------------------------------------
    // File type parameter
    //--------------------------------------------------

    std::string type = parameters_->value_string("type","unknown");

    // Error if Output::type is not defined
    if (type == "unknown") {
      ERROR1("Simulation::initialize_output_",
	     "Output:%s:type parameter is undefined",
	     file_group.c_str());
    }

    // Create output object

    Output * output = create_output_(type);

    // Error if output type was not recognized
    if (output == NULL) {
      ERROR2("Simulation::initialize_output_",
	     "Unrecognized parameter value Output:%s:type = %s",
	     file_group.c_str(),
	     type.c_str());
    }

    //--------------------------------------------------
    // File name parameter
    //--------------------------------------------------

    std::string file_name = "";
    std::vector<std::string> file_args;

    if (parameters_->type("name") == parameter_string) {

      // Case 1: string e.g. name = "filename";

      file_name = parameters_->value_string("name","");

     
    } else if (parameters_->type("name") == parameter_list) {
      // Case 2: list e.g. name = ["filename-cycle-%0d.time-%5.3f", "cycle","time"]

      int list_length = parameters_->list_length("name");

      // get file name
      if (list_length > 0) {
	file_name = parameters_->list_value_string(0,"name","");
      }

      // get file args ("cycle", "time", etc.) to schedule
      for (int index = 1; index<list_length; index++) {
	file_args.push_back(parameters_->list_value_string(index,"name",""));
      }

    } else {

      ERROR1("Simulation::initialize_output_",
	     "Bad type %d for Output 'name' parameter",parameters_->type("name"));

    }

    // Error check name specified

    ASSERT1("Simulation::initialize_output_",
	   "Output 'name' must be specified for file group %s",
	   file_group.c_str(),
	   file_name != "");

    // Set the output object file name and arguments

    output->set_filename (file_name,file_args);

    //--------------------------------------------------
    // Field_list parameter
    //--------------------------------------------------

    if (parameters_->type("field_list") != parameter_unknown) {

      // Set field list to specified field list

      if (parameters_->type("field_list") == parameter_list) {

	ItFieldList * it_field = new ItFieldList;

	int length = parameters_->list_length("field_list");
	for (int i=0; i<length; i++) {
	  std::string field_name = 
	    parameters_->list_value_string(i,"field_list","");
	  int field_index = field_descr_->field_id(field_name);
	  it_field->append(field_index);
	}
	
	output->set_it_field(it_field);

      } else {

	ERROR("Simulation::initialize_output_",
	      "Bad type for Output 'field_list' parameter");

      }

    } else {

      // field_list not defined: default to all fields
      int field_count = field_descr_->field_count();
      ItFieldRange * it_field = new ItFieldRange(field_count);

      output->set_it_field(it_field);

    }

    //--------------------------------------------------
    // Scheduling parameters
    //--------------------------------------------------

    // error check schedule parameter exists

    ASSERT("Simulation::initialize_output_",
	   "The 'schedule' parameter must be specified for all Output file groups",
	   parameters_->type("schedule") != parameter_unknown);

    // get schedule variable ("cycle" or "time")

    bool var_cycle,var_time;
    var_cycle = (strcmp(parameters_->list_value_string(0,"schedule"),"cycle")==0);
    var_time  = (strcmp(parameters_->list_value_string(0,"schedule"),"time")==0);

    // error check variable name

    ASSERT("Simulation::initialize_output_",
	   "The first 'schedule' parameter list element must be 'cycle' or 'time'",
	   var_cycle || var_time);

    // get schedule type (interval or list)

    bool type_interval,type_list;

    type_interval = (strcmp(parameters_->list_value_string(1,"schedule"),"interval")==0);
    type_list     = (strcmp(parameters_->list_value_string(1,"schedule"),"list")==0);

    // error check schedule type

    ASSERT("Simulation::initialize_output_",
	   "The second 'schedule' parameter list element "
	   "must be 'interval' or 'list'",
	   type_interval || type_list);

    int len = parameters_->list_length("schedule");

    if (var_cycle && type_interval) {

      // get cycle interval schedule

      const int max_int = std::numeric_limits<int>::max();

      // Set time limits

      int start, stop, step;

      if (len == 3) { 

	// schedule = [ "cycle", "interval",  <step> ];

	start = 0;
	stop  = max_int;
	step  = parameters_->list_value_integer(2,"schedule");

      } else if (len == 5) {

	// schedule = [ "cycle", "interval",  <start>, <stop>, <step> ]

	start = parameters_->list_value_integer(2,"schedule");
	stop  = parameters_->list_value_integer(3,"schedule");
	step  = parameters_->list_value_integer(4,"schedule");

      } else {

	// error check schedule parameter list length

	ERROR("Simulation::initialize_output_",
	      "Output 'schedule' list parameter has wrong number of elements");

      }

      // Set cycle interval output schedule 

      output->schedule()->set_cycle_interval(start,step,stop);

    } else if (var_cycle && type_list) {

      // get cycle list schedule

      std::vector<int> list;

      for (int index = 2; index < len; index++) {

	int value = parameters_->list_value_integer(index,"schedule");

	list.push_back (value);

	if (list.size() > 1) {

	  // error check monotonicity

	  ASSERT("Simulation::initialize_output_",
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

      double start, stop, step;

      if (len == 3) {

	// schedule = [ "time", "interval",  <step> ];

	start = 0;
	stop  = max_double;
	step  = parameters_->list_value_float(2,"schedule");

      } else if (len == 5) {

	// schedule = [ "time", "interval",  <start>, <stop>, <step> ]

	start = parameters_->list_value_float(2,"schedule");
	stop  = parameters_->list_value_float(3,"schedule");
	step  = parameters_->list_value_float(4,"schedule");

      } else {

	// error check schedule parameter list length

	ERROR("Simulation::initialize_output_",
	      "Output 'schedule' list parameter has wrong number of elements");

      }

      // Set time interval output schedule 

      output->schedule()->set_time_interval(start,step,stop);

    } else if (var_time && type_list) {

      // get time list schedule

      std::vector<double> list;

      for (int index = 2; index < len; index++) {

	double value = parameters_->list_value_float(index,"schedule");

	list.push_back (value);

	if (list.size() > 1) {

	  // error check monotonicity

	  ASSERT("Simulation::initialize_output_",
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

    if (output != NULL) {

      // Parse image-specific parameters

      parameter_enum type = parameter_unknown;

      // position parameter

      type = parameters_->type("position");

      if ( type != parameter_unknown) {

	// error check position type

	ASSERT1("Simulation::initialize_output_",
	       "Output %s position must be a float",
		file_group.c_str(),
		type == parameter_float);

	INCOMPLETE("Simulation::initialize_output_(): image position not implemented");

      }

      // axis parameter

      type = parameters_->type("axis");

      if (type != parameter_unknown) {

	// error check position type

	ASSERT1("Simulation::initialize_output_",
	       "Output %s axis must be a string",
	       file_group.c_str(), type == parameter_string);

	std::string axis = parameters_->value_string("axis");

	// set the output axis

	if (axis == "x") output_image->set_axis(axis_x);
	if (axis == "y") output_image->set_axis(axis_y);
	if (axis == "z") output_image->set_axis(axis_z);

	// error check axis

	ASSERT2("Simulation::initialize_output_",
	      "Output %s axis %d must be \"x\", \"y\", or \"z\"",
		file_group.c_str(), axis.c_str(),
		axis=="x" || axis=="y" || axis=="z");


      }

      // colormap parameter

      type = parameters_->type("colormap");

      if (type != parameter_unknown) {

	// error check colormap list type

	ASSERT1("Simulation::initialize_output_",
	       "Output %s colormap must be a list",
	       file_group.c_str(), type == parameter_list);

	int n = parameters_->list_length("colormap");

	// error check colormap list length

	ASSERT1("Simulation::initialize_output_",
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

	  ASSERT1("Simulation::initialize_output_",
		  "Output %s colormap list must only contain floats",
		  file_group.c_str(), 
		  ((parameters_->list_type(ir,"colormap") == parameter_float) &&
		   (parameters_->list_type(ig,"colormap") == parameter_float) &&
		   (parameters_->list_type(ib,"colormap") == parameter_float)));

	  // get next colormap r[i] g[i], b[i]
	  r[i] = parameters_->list_value_float (ir, "colormap",0.0);
	  g[i] = parameters_->list_value_float (ig, "colormap",0.0);
	  b[i] = parameters_->list_value_float (ib, "colormap",0.0);

	}

	// set the colormap

	output_image->set_colormap(n,r,g,b);

	// deallocate arrays

	delete r;
	delete g;
	delete b;

      }
    }

    // Add the initialized Output object to the Simulation's list of output objects

    output_list_.push_back(output); 

  } // (for index_file_group)

}

//----------------------------------------------------------------------

void Simulation::initialize_method_() throw()
{
  //--------------------------------------------------
  parameters_->group_set(0,"Method");
  //--------------------------------------------------

  int method_count = parameters_->list_length("sequence");

  if (method_count == 0) {
    ERROR ("Simulation::initialize_method_",
	   "List parameter 'Method sequence' must have length "
	   "greater than zero");
  }

  for (int i=0; i<method_count; i++) {

    //--------------------------------------------------
    // parameter: Method:sequence
    //--------------------------------------------------

    std::string method_name = parameters_->list_value_string(i,"sequence");

    Method * method = create_method_(method_name);
    if (method) {

      method_list_.push_back(method); 

    } else {
      ERROR1("Simulation::initialize_method_",
	     "Unknown Method %s",method_name.c_str());
    }
  }
}

//----------------------------------------------------------------------

void Simulation::initialize_parallel_() throw()
{
  //--------------------------------------------------
  // parameter: parallel::temp_update_all
  // parameter: parallel::temp_update_full
  //--------------------------------------------------

  temp_update_all_  = 
    parameters_->value_logical ("Parallel:temp_update_all",true);
  temp_update_full_ = 
    parameters_->value_logical ("Parallel:temp_update_full",true);

}
//----------------------------------------------------------------------

void Simulation::deallocate_() throw()
{
  delete factory_;       factory_     = 0;
  delete parameters_;    parameters_  = 0;
  delete performance_;   performance_ = 0;
#ifdef CONFIG_USE_CHARM
  delete monitor_;       monitor_ = 0;
  delete group_process_; group_process_ = 0;
#endif
  delete hierarchy_;     hierarchy_ = 0;
  delete field_descr_;   field_descr_ = 0;
  delete stopping_;      stopping_ = 0;
  delete timestep_;      timestep_ = 0;
  delete initial_;       initial_ = 0;
  delete boundary_;      boundary_ = 0;
  for (size_t i=0; i<output_list_.size(); i++) {
    delete output_list_[i];    output_list_[i] = 0;
  }
  for (size_t i=0; i<method_list_.size(); i++) {
    delete method_list_[i];    method_list_[i] = 0;
  }
}

//----------------------------------------------------------------------

void Simulation::run() throw()
{
  ERROR ("Simulation::run","Implictly abstract function called");
}

//----------------------------------------------------------------------

void Simulation::read() throw()
{
  ERROR ("Simulation::read","Implictly abstract function called");
}

//----------------------------------------------------------------------

void Simulation::write() const throw()
{
  ERROR ("Simulation::write","Implictly abstract function called");
}

//----------------------------------------------------------------------

const Factory & Simulation::factory() const throw()
{
  if (factory_ == NULL) factory_ = new Factory;
  return *factory_;
}

//----------------------------------------------------------------------

Stopping * Simulation::create_stopping_ (std::string name) throw ()
{
  ERROR ("Simulation::create_stopping_","Implictly abstract function called");
}

//----------------------------------------------------------------------

Timestep * Simulation::create_timestep_ (std::string name) throw ()
{ 
  ERROR ("Simulation::create_timestep_","Implictly abstract function called");
}

//----------------------------------------------------------------------

Initial * Simulation::create_initial_ (std::string name) throw ()
{ 
  ERROR ("Simulation::create_initial_","Implictly abstract function called");
}

//----------------------------------------------------------------------

Boundary * Simulation::create_boundary_ (std::string name) throw ()
{
  ERROR ("Simulation::create_boundary_","Implictly abstract function called");
}

//----------------------------------------------------------------------

Output * Simulation::create_output_ (std::string type) throw ()
{ 
  Output * output = NULL;
  if (type == "image") {
    output = new OutputImage (this);
  } else if (type == "data") {
    output = new OutputData (this);
  }
  return output;
}

//----------------------------------------------------------------------

Method * Simulation::create_method_ (std::string name) throw ()
{
  ERROR ("Simulation::create_method_","Implictly abstract function called");
}

//======================================================================


#ifdef CONFIG_USE_CHARM

Simulation::Simulation() 
{
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

Simulation::Simulation (CkMigrateMessage *m) 
{
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Simulation::refresh() throw()
{

  //--------------------------------------------------
  // Monitor
  //--------------------------------------------------

  monitor_-> print("[Simulation %d] cycle %04d time %15.12f dt %15.12g", 
		   index_,cycle_,time_,dt_);

  //--------------------------------------------------
  // Output
  //--------------------------------------------------

  for (size_t index_output=0; index_output<num_output(); index_output++) {
    if (output(index_output)->schedule()->write_this_cycle(cycle_, time_)) {
      // output_open(index_output);
    }
  }

  //--------------------------------------------------
  // Stopping
  //--------------------------------------------------

  if (stop_) {

    performance_->stop();
    performance_->print(monitor_);
    proxy_main.p_exit(CkNumPes());

  } else {

    //--------------------------------------------------
    // Boundary
    //--------------------------------------------------

    ItPatch it_patch(hierarchy_);
    Patch * patch;
    axis_enum axis = (temp_update_all_) ? axis_all : axis_x;
    while (( patch = ++it_patch )) {
      if (patch->blocks_allocated()) {
#ifdef ORIGINAL_REFRESH
	patch->block_array().p_refresh(cycle_, time_, dt_,axis);
#else
	patch->block_array().p_compute(dt_,axis);
#endif
      }
    }
  }
}

#endif

//----------------------------------------------------------------------
// NOT CHARM
//----------------------------------------------------------------------

#ifndef CONFIG_USE_CHARM

void Simulation::scheduled_output()

{
  for (size_t i=0; i<output_list_.size(); i++) {
    Output * output = output_list_[i];
    if (output->is_scheduled(cycle_,time_)) {
      output->init();
      output->open();
      output->write_hierarchy(field_descr_, hierarchy_);
      output->close();
      output->finalize();
    }
  }
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

#  include "simulation.def.h"

#endif /* CONFIG_USE_CHARM */

//======================================================================
