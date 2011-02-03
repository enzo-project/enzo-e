// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "simulation.hpp"

/// Initialize the Simulation object
Simulation::Simulation(Error      * error,
		       Monitor    * monitor,
		       Parameters * parameters)
  : dimension_(0),
    error_      (error      ? error      : new Error),
    monitor_    (monitor    ? monitor    : new Monitor),
    parameters_ (parameters ? parameters : new Parameters (monitor)),
    mesh_(0),
    data_descr_(0),
    control_(0),
    timestep_(0),
    initial_(0),
    boundary_(0),
    method_list_()
{
  // Initialize parameter defaults

  extent_[0] = 0.0;
  extent_[1] = 1.0;
  extent_[2] = 0.0;
  extent_[3] = 1.0;
  extent_[4] = 0.0;
  extent_[5] = 1.0;

}

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
{
  delete data_descr_;
  delete parameters_;
}

//----------------------------------------------------------------------

Simulation::Simulation(const Simulation & simulation) throw()
{
  INCOMPLETE_MESSAGE("Simulation::Simulation(simulation)","");
}

//----------------------------------------------------------------------


Simulation & Simulation::operator= (const Simulation & simulation) throw()
{
  INCOMPLETE_MESSAGE("Simulation::operatior = (simulation)","");
  return (*this);
}

//----------------------------------------------------------------------

void Simulation::initialize(FILE * fp) throw()
{

  // Read in parameters

  parameters_->read(fp); // MEMORY LEAK

  // Initialize parameters

  initialize_simulation_();

  // Initialize simulation components

  initialize_data_();
  initialize_mesh_();

  initialize_control_();
  initialize_timestep_();
  initialize_initial_();
  initialize_boundary_();
  initialize_method_();

}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  delete data_descr_;  data_descr_ = 0;
  delete mesh_;  mesh_ = 0;
  
}

//----------------------------------------------------------------------

int Simulation::dimension() throw()
{
  return dimension_;
}

//----------------------------------------------------------------------

void Simulation::extents (double * xmin, double *xmax,
			  double * ymin, double *ymax,
			  double * zmin, double *zmax) throw()
{
  if (xmin) *xmin = extent_[0];
  if (xmax) *xmax = extent_[1];
  if (ymin) *ymin = extent_[2];
  if (ymax) *ymax = extent_[3];
  if (zmin) *zmin = extent_[4];
  if (zmax) *zmax = extent_[5];
}

//----------------------------------------------------------------------

Error * Simulation::error() const throw()
{
  return error_;
}

//----------------------------------------------------------------------

Monitor * Simulation::monitor() const throw()
{
  return monitor_;
}

//----------------------------------------------------------------------

Mesh * Simulation::mesh() const throw()
{
  return mesh_;
}
  
//----------------------------------------------------------------------

DataDescr * Simulation::data_descr() const throw()
{
  return data_descr_;
}

//----------------------------------------------------------------------

Control * Simulation::control() const throw() 
{ return control_; }

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

int Simulation::num_method() const throw()
{ return method_list_.size(); }

//----------------------------------------------------------------------

Method * Simulation::method(int i) const throw()
{ return method_list_[i]; }

//======================================================================

void Simulation::initialize_simulation_() throw()
{

  //--------------------------------------------------
  parameters_->set_current_group("Physics");
  //--------------------------------------------------

  dimension_ = parameters_->value_integer("dimensions",0);


  parameters_->set_current_group("Domain");

  // get extent_length and check consistency

  int extent_length = parameters_->list_length("extent");

  bool valid_extent_length = ((extent_length == 2) ||
			      (extent_length == 4) ||
			      (extent_length == 6));
  ASSERT ("Simulation::initialize_simulation_",
	  "Parameter Domain:extent list must have length 2, 4, or 6",
	  valid_extent_length);

  ASSERT ("Simulation::initialize_simulation_",
	  "Parameter mismatch between Physics:dimensions and Domain:extent",
	  (dimension_ == extent_length / 2));

  for (int i=0; i<extent_length; i++) {

    double extent_default = (i % 2 == 0) ? 0.0 : 1.0;

    extent_[i] = parameters_->list_value_scalar(i, "extent", extent_default);

    if (i % 2 == 1) {
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,
	       "Parameter Domain:extent[%g] not lower than Domain:extent[%g]",
	       extent_[i-1],extent_[i]);
      ASSERT ("Simulation::Simulation", error_message, (extent_[i-1] < extent_[i]));
    }
  }

}

//----------------------------------------------------------------------

void Simulation::initialize_data_() throw()
{

  data_descr_ = new DataDescr;

  FieldDescr * field_descr = data_descr_->field_descr();

  parameters_->set_current_group("Field");

  // Add data fields

  int i;
  for (i=0; i<parameters_->list_length("fields"); i++) {
    field_descr->insert_field
      (parameters_->list_value_string(i, "fields"));
  }

  // Define default ghost zone depth for all fields, default value of 1

  int gx = parameters_->list_value_integer(0,"ghosts",1);
  int gy = parameters_->list_value_integer(1,"ghosts",1);
  int gz = parameters_->list_value_integer(2,"ghosts",1);

  if (dimension_ < 2) gy = 0;
  if (dimension_ < 3) gz = 0;

  for (i=0; i<field_descr->field_count(); i++) {
    field_descr->set_ghosts(i,gx,gy,gz);
  }

}

//----------------------------------------------------------------------

void Simulation::initialize_mesh_() throw()
{

  ASSERT("Simulation::initialize_mesh_",
	 "data must be initialized before mesh",
	 data_descr_ != NULL);

  mesh_ = new Mesh(data_descr_);

  parameters_->set_current_group ("Mesh");

  // Block sizes

  int block_size[3];

  block_size[0] = parameters_->list_value_integer(0,"block_size",4);
  block_size[1] = parameters_->list_value_integer(1,"block_size",4);
  block_size[2] = parameters_->list_value_integer(2,"block_size",4);

  // assume constant block size
  mesh_->set_min_block_size (block_size[0]);
  mesh_->set_max_block_size (block_size[0]);
  
  // Patch sizes

  int patch_size[3];

  patch_size[0] = parameters_->list_value_integer(0,"patch_size",128);
  patch_size[1] = parameters_->list_value_integer(1,"patch_size",128);
  patch_size[2] = parameters_->list_value_integer(1,"patch_size",128);

  // minimum patch size is one block
  mesh_->set_min_patch_size (block_size[0]);
  mesh_->set_max_patch_size (patch_size[0]);

  // Mesh size

  int root_size[3];

  root_size[0] = parameters_->list_value_integer(0,"root_size",1);
  root_size[1] = parameters_->list_value_integer(1,"root_size",1);
  root_size[2] = parameters_->list_value_integer(2,"root_size",1);

  mesh_->set_root_size(root_size[0],root_size[1],root_size[2]);

  parameters_->set_current_group("Domain");

  double domain_lower[3];

  domain_lower[0] = parameters_->list_value_scalar(0,"extent",0.0);
  domain_lower[1] = parameters_->list_value_scalar(2,"extent",0.0);
  domain_lower[2] = parameters_->list_value_scalar(4,"extent",0.0);

  mesh_->set_lower(domain_lower[0],
		   domain_lower[1],
		   domain_lower[2]);

  double domain_upper[3];

  domain_upper[0] = parameters_->list_value_scalar(1,"extent",1.0);
  domain_upper[1] = parameters_->list_value_scalar(3,"extent",1.0);
  domain_upper[2] = parameters_->list_value_scalar(5,"extent",1.0);

  mesh_->set_upper(domain_upper[0],
		   domain_upper[1],
		   domain_upper[2]);

  // Parallel layout of the root mesh patch

  
  Layout * layout = new Layout;

  parameters_->set_current_group("Mesh");

  int process_first = parameters_->value_integer("root_process_first",0);
  int process_count = parameters_->value_integer("root_process_count",1);

  layout->set_process_range(process_first, process_count);

  int block_count[3];

  block_count[0] = parameters_->list_value_integer(0,"root_blocks",1);
  block_count[1] = parameters_->list_value_integer(1,"root_blocks",1);
  block_count[2] = parameters_->list_value_integer(2,"root_blocks",1);

  layout->set_block_count(block_count[0],block_count[1],block_count[2]);
    

  // Create the root patch

  Patch * root_patch = new Patch;

  root_patch->set_data_descr(data_descr_);
  root_patch->set_size(root_size[0],root_size[1],root_size[2]);
  root_patch->set_layout(layout);
  root_patch->set_extents(domain_lower[0],domain_upper[0],
			  domain_lower[1],domain_upper[1],
			  domain_lower[2],domain_upper[2]);

  mesh_->set_root_patch(root_patch);

  // Mesh AMR settings

  mesh_->set_refine    (parameters_->value_integer("refine",    2));
  mesh_->set_max_level (parameters_->value_integer("max_level", 0));
  mesh_->set_balanced  (parameters_->value_logical("balanced",  true));
  mesh_->set_backfill  (parameters_->value_logical("backfill",  true));
  mesh_->set_coalesce  (parameters_->value_logical("coalesce",  true));

}

//----------------------------------------------------------------------

void Simulation::initialize_control_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Control");
  //--------------------------------------------------

  std::string name = parameters_->value_string("name","");
  control_ = create_control_(name);
}

//----------------------------------------------------------------------

void Simulation::initialize_timestep_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Timestep");
  //--------------------------------------------------

  std::string name = parameters_->value_string("name","");
  timestep_ = create_timestep_(name);
}

//----------------------------------------------------------------------

void Simulation::initialize_initial_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Initial");
  //--------------------------------------------------

  std::string name = parameters_->value_string("name","");
  initial_ = create_initial_(name);
}

//----------------------------------------------------------------------

void Simulation::initialize_boundary_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Boundary");
  //--------------------------------------------------

  std::string name = parameters_->value_string("name","");
  boundary_ = create_boundary_(name);
}


//----------------------------------------------------------------------

void Simulation::initialize_method_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Method");
  //--------------------------------------------------

  int method_count = parameters_->list_length("sequence");

  if (method_count == 0) {
    ERROR_MESSAGE ("Simulation::initialize_method_",
		   "List parameter 'Method sequence' must have length greater than zero");
  }

  for (int i=0; i<method_count; i++) {

    std::string method_name = parameters_->list_value_string(i,"sequence");

    Method * method = create_method_(method_name);
    if (method) {

      method_list_.push_back(method); 
      method->initialize(data_descr_);
    } else {
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"Unknown Method %s",method_name.c_str());
      ERROR_MESSAGE ("Simulation::initialize_method_",
		     error_message);
    }
  }


  INCOMPLETE_MESSAGE("Simulation::initialize_method_","");
}

