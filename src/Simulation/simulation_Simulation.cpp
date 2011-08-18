// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class
/// @todo      Move dt / cycle update to separate function from p_output()

#include "cello.hpp"

class Factory;
#include "main.hpp"

#include "simulation.hpp"
#include "simulation_charm.hpp"

Simulation::Simulation
(
 const char *   parameter_file,
#ifdef CONFIG_USE_CHARM
 int n,
 const Factory & factory,
 CProxy_BlockReduce proxy_block_reduce,
#else
 const Factory & factory,
 GroupProcess * group_process,
#endif
 int            index
 ) throw()
/// Initialize the Simulation object
  : parameters_(0),
    factory_(&factory),
#ifndef CONFIG_USE_CHARM
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
    mesh_(0),
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
  parameters_  = new Parameters(parameter_file,monitor_);
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
  initialize_mesh_();
  initialize_stopping_();
  initialize_timestep_();
  initialize_initial_();
  initialize_boundary_();
  initialize_output_();
  initialize_method_();
  initialize_parallel_();

  
  char parameter_file[40];
  snprintf (parameter_file,40,"cello-%d.%g.out",cycle_,time_);
  parameters_->write(parameter_file);
}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

int Simulation::dimension() const throw()
{
  return dimension_;
}

//----------------------------------------------------------------------

Mesh * Simulation::mesh() const throw()
{
  return mesh_;
}
  
//----------------------------------------------------------------------

FieldDescr * Simulation::field_descr() const throw()
{
  return field_descr_;
}

//----------------------------------------------------------------------

Performance * Simulation::performance() const throw()
{
  return performance_;
}

//----------------------------------------------------------------------

Monitor * Simulation::monitor() const throw()
{
  return monitor_;
}

//----------------------------------------------------------------------

Stopping * Simulation::stopping() const throw() 
{ return stopping_; }

//----------------------------------------------------------------------

Timestep * Simulation::timestep() const throw() 
{ return timestep_; }

//----------------------------------------------------------------------

Initial * Simulation::initial() const throw() 
{ return initial_; }

//----------------------------------------------------------------------

Boundary * Simulation::boundary() const throw() 
{ return boundary_; }

//----------------------------------------------------------------------

size_t Simulation::num_output() const throw()
{ return output_list_.size(); }

//----------------------------------------------------------------------

Output * Simulation::output(int i) const throw()
{ return output_list_[i]; }

//----------------------------------------------------------------------

size_t Simulation::num_method() const throw()
{ return method_list_.size(); }

//----------------------------------------------------------------------

Method * Simulation::method(int i) const throw()
{ return method_list_[i]; }

//----------------------------------------------------------------------

const Factory * Simulation::factory () const throw()
{ return factory_; }

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

void Simulation::initialize_mesh_() throw()
{

  ASSERT("Simulation::initialize_mesh_",
	 "data must be initialized before mesh",
	 field_descr_ != NULL);

  //----------------------------------------------------------------------
  // Create and initialize Mesh
  //----------------------------------------------------------------------

  mesh_ = factory()->create_mesh ();

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

  mesh_->set_lower(lower[0], lower[1], lower[2]);
  mesh_->set_upper(upper[0], upper[1], upper[2]);

  //----------------------------------------------------------------------
  // Create and initialize root Patch in Mesh
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

  mesh_->create_root_patch
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

  int num_groups = parameters_->list_length("groups");

  for (int index_group=0; index_group < num_groups; index_group++) {

    //--------------------------------------------------
    parameters_->group_set(0,"Output");
    //--------------------------------------------------

    std::string group = parameters_->list_value_string
      (index_group,"groups","unknown");

    //--------------------------------------------------
    parameters_->group_set(1,group);
    //--------------------------------------------------
    // parameter: Output:<group>:type
    // parameter: Output:<group>:file_name
    // parameter: Output:<group>:field_list
    // parameter: Output:<group>:cycle_interval
    // parameter: Output:<group>:cycle_list
    // parameter: Output:<group>:time_interval
    // parameter: Output:<group>:time_list
    //--------------------------------------------------

    std::string type = parameters_->value_string("type","unknown");

    // Error if Output::type is not defined
    if (type == "unknown") {
      char buffer[ERROR_LENGTH];
      sprintf (buffer,"Output:%s:type parameter is undefined",group.c_str());
      ERROR("Simulation::initialize_output_",buffer);
    }

    Output * output = create_output_(type);

    // Error if output type was not recognized
    if (output == NULL) {
      char buffer[ERROR_LENGTH];
      sprintf (buffer,"Unrecognized parameter value Output:%s:type = %s",
     	       group.c_str(),type.c_str());
      ERROR("Simulation::initialize_output_",buffer);
    }

    // ASSUMES GROUP AND SUBGROUP ARE SET BY CALLER

    ASSERT("Simulation::initialize_output_",
	   "Bad type for Output 'file_name' parameter",
	   parameters_->type("file_name") == parameter_string);

    // Set the file name

    std::string file_name = parameters_->value_string("file_name","unknown");

    // Error if file_name is unspecified
    ASSERT("Simulation::initialize_output_",
	   "Output 'file_name' must be specified",
	   file_name != "unknown");

    output->set_file_name (file_name);

    std::vector<int> field_list;

    if (parameters_->type("field_list") != parameter_unknown) {
      // Set field list to specified field list
      if (parameters_->type("field_list") == parameter_list) {
	int length = parameters_->list_length("field_list");
	for (int i=0; i<length; i++) {
	  int field_index = parameters_->list_value_integer(i,"field_list",0);
	  field_list.push_back(field_index);
	}
	output->set_field_list(field_list);
      } else {
	ERROR("Simulation::initialize_output_",
	      "Bad type for Output 'field_list' parameter");
      }
    } else {
      // Set field list to default of all fields
      int field_count = field_descr_->field_count();
      for (int i=0;  i<field_count; i++) {
	field_list.push_back(i);
      }
      output->set_field_list(field_list);
    }

    //--------------------------------------------------
    // Determine scheduling
    //--------------------------------------------------

    bool cycle_interval,cycle_list,time_interval,time_list;

    cycle_interval = (parameters_->type("cycle_interval") != parameter_unknown);
    cycle_list     = (parameters_->type("cycle_list")     != parameter_unknown);
    time_interval  = (parameters_->type("time_interval")  != parameter_unknown);
    time_list      = (parameters_->type("time_list")      != parameter_unknown);

    ASSERT("Simulation::initialize_output_",
	   "exactly one of [cycle|time]_[interval|list] must be defined for each Output group",
	   (cycle_interval? 1 : 0 + 
	    cycle_list?     1 : 0 + 
	    time_interval?  1 : 0 + 
	    time_list?      1 : 0) == 1);

    if (cycle_interval) {

      const int max_int = std::numeric_limits<int>::max();
      int cycle_start = 0;
      int cycle_step = 0;
      int cycle_stop = max_int;
      if (parameters_->type("cycle_interval") == parameter_integer) {
	cycle_step = parameters_->value_integer("cycle_interval",1);
      } else if (parameters_->type("cycle_interval") == parameter_list) {
	if (parameters_->list_length("cycle_interval") != 3) {
	  ERROR("Simulation::initialize_output_",
		"Output cycle_interval parameter must be of the form "
		"[cycle_start, cycle_step, cycle_stop");
	}
	cycle_start = 
	  parameters_->list_value_integer (0,"cycle_interval",0);
	cycle_step  = 
	  parameters_->list_value_integer (1,"cycle_interval",1);
	cycle_stop  = 
	  parameters_->list_value_integer (2,"cycle_interval",max_int);
      } else {
	ERROR("Simulation::initialize_output_",
	      "Output cycle_interval is of the wrong type");
      }
      output->set_cycle_interval(cycle_start,cycle_step,cycle_stop);

    } else if (cycle_list) {

      std::vector<int> list;

      if (parameters_->type("cycle_list") == parameter_integer) {
	list.push_back (parameters_->value_integer("cycle_list",0));
      } else if (parameters_->type("cycle_list") == parameter_list) {
	for (int i=0; i<parameters_->list_length("cycle_list"); i++) {
	  int value = parameters_->list_value_integer(i,"cycle_list",0);
	  list.push_back (value);
	  // check monotonicity
	  if (i > 0) {
	    int value_prev = parameters_->list_value_integer(i-1,"cycle_list",0);
	    ASSERT("Simulation::initialize_output_",
		   "Output cycle_list must be monotonically increasing",
		   value_prev < value);
	  }
	}
      }

      output->set_cycle_list(list);

    } else if (time_interval) {

      const double max_double = std::numeric_limits<double>::max();
      double time_start = 0;
      double time_step = 0;
      double time_stop = max_double;
      if (parameters_->type("time_interval") == parameter_float) {
	time_step = parameters_->value_float("time_interval",1);
      } else if (parameters_->type("time_interval") == parameter_list) {
	if (parameters_->list_length("time_interval") != 3) {
	  ERROR("Simulation::initialize_output_",
		"Output time_interval parameter must be of the form "
		"[time_start, time_step, time_stop");
	}
	time_start = 
	  parameters_->list_value_float (0,"time_interval",0);
	time_step  = 
	  parameters_->list_value_float (1,"time_interval",1);
	time_stop  = 
	  parameters_->list_value_float (2,"time_interval",max_double);
      } else {
	ERROR("Simulation::initialize_output_",
	      "Output time_interval is of the wrong type");
      }
      output->set_time_interval(time_start,time_step,time_stop);

    } else if (time_list) {

      std::vector<double> list;

      if (parameters_->type("time_list") == parameter_float) {
	list.push_back (parameters_->value_float("time_list",0));
      } else if (parameters_->type("time_list") == parameter_list) {
	for (int i=0; i<parameters_->list_length("time_list"); i++) {
	  double value = parameters_->list_value_float(i,"time_list",0.0);
	  list.push_back (value);
	  // check monotonicity
	  if (i > 0) {
	    double value_prev = parameters_->list_value_float(i-1,"time_list",0);
	    ASSERT("Simulation::initialize_output_",
		   "Output time_list must be monotonically increasing",
		   value_prev < value);
	  }
	}
      }

      output->set_time_list(list);

    }

    output_list_.push_back(output); 

  } // (for index_group)

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
	   "List parameter 'Method sequence' must have length greater than zero");
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
      char error_message[ERROR_LENGTH];
      sprintf (error_message,"Unknown Method %s",method_name.c_str());
      ERROR ("Simulation::initialize_method_",
		     error_message);
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
  delete performance_;
  performance_ = 0;
#ifdef CONFIG_USE_CHARM
  delete monitor_;
  monitor_ = 0;
#endif
  delete mesh_;
  mesh_ = 0;
  delete field_descr_;
  field_descr_ = 0;
  delete stopping_;
  stopping_ = 0;
  delete timestep_;
  timestep_ = 0;
  delete initial_;
  initial_ = 0;
  delete boundary_;
  boundary_ = 0;
  delete factory_;
  factory_ = 0;
  for (size_t i=0; i<method_list_.size(); i++) {
    delete method_list_[i];
    method_list_[i] = 0;
  }
}

void Simulation::run() throw()
{ ERROR ("Simulation::run","Implictly abstract function called"); }

void Simulation::read() throw()
{ ERROR ("Simulation::read","Implictly abstract function called"); }

void Simulation::write() const throw()
{ ERROR ("Simulation::write","Implictly abstract function called"); }

Stopping * Simulation::create_stopping_ (std::string name) throw ()
{ ERROR ("Simulation::create_stopping_","Implictly abstract function called"); }

Timestep * Simulation::create_timestep_ (std::string name) throw ()
{ ERROR ("Simulation::create_timestep_","Implictly abstract function called"); }

Initial * Simulation::create_initial_ (std::string name) throw ()
{ ERROR ("Simulation::create_initial_","Implictly abstract function called"); }

Boundary * Simulation::create_boundary_ (std::string name) throw ()
{ ERROR ("Simulation::create_boundary_","Implictly abstract function called"); }

Output * Simulation::create_output_ (std::string name) throw ()
{ ERROR ("Simulation::create_output_","Implictly abstract function called"); }

Method * Simulation::create_method_ (std::string name) throw ()
{ ERROR ("Simulation::create_method_","Implictly abstract function called"); }

//======================================================================

#ifdef CONFIG_USE_CHARM

void Simulation::p_output 
( int cycle, double time, double dt, bool stop ) throw()
{

  TRACE("Simulation::p_output");

  // Update Simulation cycle and time from reduction to main
  
  cycle_ = cycle;
  time_  = time;
  dt_    = dt;
  stop_  = stop;

  // reset output "loop" over output objects
  output_first();

  // process first output object, which continues with refresh() if done
  output_next();
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Simulation::output_first() throw()
{
  TRACE("Simulation::output_first()");
  index_output_ = 0;
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Simulation::output_next() throw()
{

  TRACE("Simulation::output_next()");

  // find next output

  while (index_output_ < num_output() && 
	 ! output(index_output_)->write_this_cycle(cycle_, time_))
    ++index_output_;

  // output if any scheduled, else proceed with refresh

  if (index_output_ < num_output()) {

    // Open the file(s)
    output(index_output_)->open(mesh_,cycle_,time_);

    // Call blocks to contribute their data
    ItPatch it_patch(mesh_);
    Patch * patch;
    while (( patch = ++it_patch )) {
      if (patch->blocks_allocated()) {
	patch->block_array().p_output (index_output_);
      }
    }

  } else {

    refresh();

  }
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Simulation::p_output_reduce() throw()
{
  TRACE("Simulation::p_output_reduce");

  int ip       = CkMyPe();
  int ip_write = ip - (ip % output(index_output_)->process_write());

  // Even self calls this to avoid hanging if case np == 1
  char buffer[20];
  sprintf(buffer,"%02d > %02d send",ip,ip_write);
  if (ip != ip_write) {
    PARALLEL_PRINTF("%d -> %d calling p_output_write()\n",ip,ip_write);
    proxy_simulation[ip_write].p_output_write (strlen(buffer),buffer);
    output_next();
  } else {
    PARALLEL_PRINTF("%d -> %d calling p_output_write()\n",ip,ip_write);
    proxy_simulation[ip].p_output_write(0,0);
  }

}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Simulation::p_output_write (int n, char * buffer) throw()
{
  TRACE("Simulation::p_output_write");

  Output * out = output(index_output_);
  int ip       = CkMyPe();
  int ip_write = ip - (ip % out->process_write());
  PARALLEL_PRINTF ("%d %d  %d  %d\n",ip,ip_write,CkMyPe(),out->process_write());

  int count = out->counter();

  if (count == 0) {
    PARALLEL_PRINTF("Initialize writer\n");
  }
  if (n == 0) {
    PARALLEL_PRINTF ("Process reduce this %d\n",ip);
  } else {
    PARALLEL_PRINTF ("Process reduce that\n");
  }
  
  if (count == out->process_write()) {

    PARALLEL_PRINTF ("File write / close / next\n");

    // write
    // close
    out->close();

    // next


    output_next();
  }


}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void Simulation::refresh() throw()
{

  TRACE("Simulation::refresh");

  //--------------------------------------------------
  // Monitor
  //--------------------------------------------------

  monitor_-> print("[Simulation %d] cycle %04d time %15.12f dt %15.12g", 
		   index_,cycle_,time_,dt_);

  //--------------------------------------------------
  // Output
  //--------------------------------------------------

  for (size_t index_output=0; index_output<num_output(); index_output++) {
    if (output(index_output)->write_this_cycle(cycle_, time_)) {
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

    ItPatch it_patch(mesh_);
    Patch * patch;
    axis_enum axis = (temp_update_all_) ? axis_all : axis_x;
    while (( patch = ++it_patch )) {
      if (patch->blocks_allocated()) {
#ifdef ORIGINAL_REFRESH
	patch->block_array().p_refresh(dt_,axis);
#else
	patch->block_array().p_compute(dt_,axis);
#endif
      }
    }
  }
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

#  include "simulation.def.h"

#endif /* CONFIG_USE_CHARM */

//======================================================================
