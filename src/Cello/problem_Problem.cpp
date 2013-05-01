// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Problem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-03
/// @brief    Implementation of the Problem container class

#include "problem.hpp"

//----------------------------------------------------------------------

Problem::Problem() throw()
  : boundary_(0),
    num_initial_(0),
    num_refine_(0),
    stopping_(0),
    timestep_(0),
    num_method_(0),
    num_output_(0),
    index_output_(0)
{
}

//----------------------------------------------------------------------

Problem::~Problem() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Problem::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  PUP::able::pup(p);

  bool up = p.isUnpacking();

  p | boundary_; // PUP::able

  p | num_initial_;
  TRACE1 ("num_initial_ = %d",num_initial_);
  if (up) initial_list_.resize(num_initial_);
  for (int i=0; i<num_initial_; i++) {
    p | initial_list_[i]; // PUP::able
  }

  p | num_refine_;
  TRACE1 ("num_refine_ = %d",num_refine_);
  if (up) refine_list_.resize(num_refine_);
  for (int i=0; i<num_refine_; i++) {
    p | refine_list_[i]; // PUP::able
  }

  if (up) stopping_ = new Stopping;
  p | *stopping_;

  p | timestep_; // PUP::able

  p | num_method_;
  TRACE1 ("num_method_ = %d",num_method_);
  if (up) method_list_.resize(num_method_);
  for (int i=0; i<num_method_; i++) {
    p | method_list_[i]; // PUP::able
  }

  p | num_output_;
  TRACE1 ("num_output_ = %d",num_output_);
  if (up) output_list_.resize(num_output_);
  for (int i=0; i<num_output_; i++) {
    p | output_list_[i]; // PUP::able
  }

  p | index_output_;

}

#endif

//----------------------------------------------------------------------

void Problem::initialize_boundary(Config * config) throw()
{

  std::string type = config->boundary_type;

  boundary_ = create_boundary_(type,config);

  ASSERT1("Problem::initialize_boundary",
	  "Boundary type %s not recognized",
	  type.c_str(),
	  boundary_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_initial(Config * config,
				 Parameters * parameters,
				 const GroupProcess * group_process) throw()
{
  Initial * initial = create_initial_
    (config->initial_type,parameters,config,group_process);

  initial_list_.push_back( initial );

  ++ num_initial_;

  ASSERT1("Problem::initialize_initial",
	  "Initial type %s not recognized",
	  config->initial_type.c_str(),
	  initial != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_refine(Config * config,
				FieldDescr * field_descr) throw()
{
  for (size_t i=0; i<config->mesh_adapt_type.size(); i++) {

    std::string name = config->mesh_adapt_type[i];

    Refine * refine = create_refine_ (name,config,field_descr,i);

    if (refine) {
      refine_list_.push_back( refine );
      ++ num_refine_;
    } else {
      ERROR1("Problem::initialize_refine",
	     "Unknown Refine %s",name.c_str());
    }
  }
}

//----------------------------------------------------------------------

void Problem::initialize_stopping(Config * config) throw()
{
  stopping_ = create_stopping_("default",config);

  ASSERT("Problem::initialize_stopping",
	  "Stopping object not successfully created",
	  stopping_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_timestep(Config * config) throw()
{
  //--------------------------------------------------
  // parameter: Timestep : type
  //--------------------------------------------------

  timestep_ = create_timestep_(config->timestep_type,config);

  ASSERT1("Problem::initialize_timestep",
	  "Timestep type %s not recognized",
	  config->timestep_type.c_str(),
	  timestep_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_output
(Config * config,
 FieldDescr * field_descr,
 const GroupProcess * group_process,
 const Factory * factory) throw()
{

  for (int index=0; index < config->num_file_groups; index++) {

    std::string file_group = config->output_file_groups [index];
    std::string type       = config->output_type[index];

    Output * output = create_output_ (type,index, config,group_process,factory);

    if (output == NULL) {
      ERROR2("Problem::initialize_output",
	     "Unknown parameter type Output:%s:type = %s",
	     file_group.c_str(),type.c_str());
    }

    if (config->output_name[index].size() > 0) {
      std::string file_name = config->output_name[index][0];

      std::vector<std::string> file_args;

      for (size_t i=1; i<config->output_name[index].size(); i++) {
	file_args.push_back(config->output_name[index][i]);
      }

      output->set_filename (file_name,file_args);
    }


    if (config->output_field_list[index].size() > 0) {

	ItFieldList * it_field = new ItFieldList;
	int length = config->output_field_list[index].size();
	for (int i=0; i<length; i++) {
	  std::string field_name = config->output_field_list[index][i];
	  int field_index = field_descr->field_id(field_name);
	  it_field->append(field_index);
	}
	
	output->set_it_field(it_field);

    } else {

      int field_count = field_descr->field_count();
      ItFieldRange * it_field = new ItFieldRange(field_count);

      output->set_it_field(it_field);

    }


    //--------------------------------------------------
    // Scheduling parameters
    //--------------------------------------------------

    bool var_cycle,var_time;

    var_cycle = (config->output_schedule_var[index] == "cycle");
    var_time  = (config->output_schedule_var[index] == "time");

    bool type_interval,type_list;

    type_interval = config->output_schedule_type[index] == "interval";
    type_list     = config->output_schedule_type[index] == "list";

    
    if (type_interval) {

      if (var_cycle) {

	int start = config->output_schedule_start[index];
	int stop  = config->output_schedule_stop[index];
	int step  = config->output_schedule_step[index];

	output->schedule()->set_cycle_interval(start,step,stop);

      } else if (var_time) {

	double start = config->output_schedule_start[index];
	double stop  = config->output_schedule_stop[index];
	double step  = config->output_schedule_step[index];

	output->schedule()->set_time_interval(start,step,stop);
      }

    } else if (type_list) {

      if (var_cycle) {

	int size = config->output_schedule_list[index].size();
	std::vector<int> list;
	list.resize(size);
	for (int i=0; i<size; i++) {
	  list[i] = (int) config->output_schedule_list[index][i];
	}
	output->schedule()->set_cycle_list(list);

      } else if (var_time) {

	output->schedule()->set_time_list(config->output_schedule_list[index]);

      }

    }

    //--------------------------------------------------
    // Image parameters
    //--------------------------------------------------

    OutputImage * output_image = dynamic_cast<OutputImage *> (output);

    if (output != NULL) {

      // AXIS

      std::string axis = config->output_image_axis[index];

      if (axis == "x") output_image->set_axis(axis_x);
      if (axis == "y") output_image->set_axis(axis_y);
      if (axis == "z") output_image->set_axis(axis_z);

      // COLORMAP

      int n = config->output_image_colormap[index].size() / 3;

      if (n > 0) {
	double * r = new double [n];
	double * g = new double [n];
	double * b = new double [n];

	for (int i=0; i<n; i++) {

	  int ir=3*i+0;
	  int ig=3*i+1;
	  int ib=3*i+2;

	  r[i] = config->output_image_colormap[index][ir];
	  g[i] = config->output_image_colormap[index][ig];
	  b[i] = config->output_image_colormap[index][ib];

	}

	output_image->set_colormap(n,r,g,b);

	delete r;
	delete g;
	delete b;

      }

      // R G B A

      // COLORMAP_ALPHA

      n = config->output_image_colormap_alpha[index].size() / 4;

      if (n > 0) {

	double * r = new double [n];
	double * g = new double [n];
	double * b = new double [n];
	double * a = new double [n];

	for (int i=0; i<n; i++) {

	  int ir=4*i+0;
	  int ig=4*i+1;
	  int ib=4*i+2;
	  int ia=4*i+3;

	  r[i] = config->output_image_colormap_alpha[index][ir];
	  g[i] = config->output_image_colormap_alpha[index][ig];
	  b[i] = config->output_image_colormap_alpha[index][ib];
	  a[i] = config->output_image_colormap_alpha[index][ia];

	}

	output_image->set_colormap(n,r,g,b,a);

	delete r;
	delete g;
	delete b;
	delete a;
      }
    }

    // Add the initialized Output object to the Simulation's list of
    // output objects

    output_list_.push_back(output); 
    ++ num_output_;

  } // (for index)

}

//----------------------------------------------------------------------

void Problem::initialize_method(Config * config) throw()
{

  for (size_t i=0; i<config->method_sequence.size(); i++) {

    std::string name = config->method_sequence[i];

    Method * method = create_method_(name);

    if (method) {

      method_list_.push_back(method); 
      ++ num_method_;

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
  num_initial_ = 0;
  for (size_t i=0; i<initial_list_.size(); i++) {
    delete initial_list_[i];    initial_list_[i] = 0;
  }
  num_refine_ = 0;
  for (size_t i=0; i<refine_list_.size(); i++) {
    delete refine_list_[i];    refine_list_[i] = 0;
  }
  delete stopping_;      stopping_ = 0;
  delete timestep_;      timestep_ = 0;
  num_output_ = 0;
  for (size_t i=0; i<output_list_.size(); i++) {
    delete output_list_[i];    output_list_[i] = 0;
  }
  num_method_ = 0;
  for (size_t i=0; i<method_list_.size(); i++) {
    delete method_list_[i];    method_list_[i] = 0;
  }
}

//----------------------------------------------------------------------

Boundary * Problem::create_boundary_
(
 std::string  name,
 Config * config
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
 Config * config,
 const GroupProcess * group_process
 ) throw ()
{ 
  //--------------------------------------------------
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

  if (type == "file" || type == "restart") {
    return new InitialFile   (parameters,group_process,config->initial_cycle,config->initial_time);;
  } else if (type == "default") {
    return new InitialDefault(parameters,config->initial_cycle,config->initial_time);
  }
  return NULL;
}

//----------------------------------------------------------------------

Refine * Problem::create_refine_
(
 std::string  type,
 Config * config,
 FieldDescr * field_descr,
 int index
 ) throw ()
{ 
  TRACE3("mesh_root_size = %d %d %d",
	 config->mesh_root_size[0],
	 config->mesh_root_size[1],
	 config->mesh_root_size[2]);

  if (type == "slope") {
    return new RefineSlope (field_descr,
			    config->mesh_adapt_slope_min_refine,
			    config->mesh_adapt_slope_max_coarsen,
			    config->mesh_adapt_fields);
  } else if (type == "mass") {
    double root_cell_volume = 1.0;
    for (int i=0; i<config->mesh_root_rank; i++) {
      root_cell_volume *= 
	(config->domain_upper[i] - config->domain_lower[i])
	/ (config->mesh_root_size[i]);
    }

    return new RefineMass (config->mesh_adapt_mass_min,
			   config->mesh_adapt_mass_level_exponent,
			   config->mesh_adapt_mass_min_overdensity,
			   root_cell_volume);
  }
  return NULL;
}

//----------------------------------------------------------------------

Stopping * Problem::create_stopping_ 
(
 std::string  type,
 Config * config
 ) throw ()
/// @param type   Type of the stopping method to create (ignored)
/// @param stop_cycle  Stopping cycle
/// @param stop_time  Stopping time
{
  return new Stopping(config->stopping_cycle,
		      config->stopping_time);
}

//----------------------------------------------------------------------

Timestep * Problem::create_timestep_ 
(
 std::string  type,
 Config * config
 ) throw ()
{ 
  // No default timestep
  return NULL;
}

//----------------------------------------------------------------------

Method * Problem::create_method_ ( std::string  name ) throw ()
{
  TRACE1("Problem::create_method %s",name.c_str());
  // No default method
  return NULL;
}


//----------------------------------------------------------------------

Output * Problem::create_output_ 
(
 std::string    name,
 int index,
 Config *  config,
 const GroupProcess * group_process,
 const Factory * factory
 ) throw ()
/// @param name          Name of Output object to create
/// @param group_process Image output needs group process size
{ 

  Output * output = NULL;

  if (name == "image") {

    //--------------------------------------------------
    // parameter: Mesh : root_size
    // parameter: Mesh : root_blocks
    //--------------------------------------------------

    int nx,ny;

    nx = config->mesh_root_size[0];
    ny = config->mesh_root_size[1];

    // NOTE: assumes cube for non-z axis images

    std::string image_type = config->output_image_type[index];

    output = new OutputImage (index,factory,group_process->size(),nx,ny,image_type);

  } else if (name == "data") {

    output = new OutputData (index,factory,config);

  } else if (name == "restart") {

    output = new OutputRestart (index,factory,config,group_process->size());

  }

  return output;

}

