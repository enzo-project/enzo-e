// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "simulation.hpp"
#include "mesh.hpp"
#include "method.hpp" 

/// Initialize the Simulation object
Simulation::Simulation(Global * global)
  : dimension_(0),
    global_(global),
    //
    mesh_(0),
    data_descr_(0),
    initialize_list_(),
    method_list_(),
    timestep_(0),
    control_(0)
{
  // Initialize parameter defaults

  extent_[0] = 0.0;
  extent_[1] = 1.0;
  extent_[2] = 0.0;
  extent_[3] = 1.0;
  extent_[4] = 0.0;
  extent_[5] = 1.0;

  // Initialize simulation component defaults

  data_descr_ = new DataDescr;

  mesh_ = new Mesh(data_descr_);

}

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
{
  delete data_descr_;
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

void Simulation::initialize(std::string parameter_file) throw()
{
  Parameters * parameters = global_->parameters();

  // Read in parameters

  FILE * file_pointer;
  file_pointer = fopen (parameter_file.c_str(),"r");
  parameters->read(file_pointer); // MEMORY LEAK
  fclose(file_pointer);

  // Initialize parameters

  initialize_simulation_(parameters);

  // Initialize simulation components

  initialize_mesh_(parameters);
  initialize_data_(parameters);

  initialize_initial_(parameters);
  initialize_method_(parameters);
  initialize_control_(parameters);
  initialize_timestep_(parameters);

}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  INCOMPLETE_MESSAGE("Simulation::finalize","");
}

//----------------------------------------------------------------------

void Simulation::run() throw()
{
  INCOMPLETE_MESSAGE("Simulation::run","");
}

//----------------------------------------------------------------------

void Simulation::read() throw()
{
  INCOMPLETE_MESSAGE("Simulation::read","");
}

//----------------------------------------------------------------------

void Simulation::write() throw()
{
  INCOMPLETE_MESSAGE("Simulation::write","");
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

Global * Simulation::global() const throw()
{
  return global_;
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

int Simulation::num_initial() const throw()
{
  return initialize_list_.size();
}

//----------------------------------------------------------------------

MethodInitial * Simulation::initial(int i) const throw()
{
  return initialize_list_[i];
}

//----------------------------------------------------------------------

int Simulation::num_method() const throw()
{
  return method_list_.size();
}

//----------------------------------------------------------------------

MethodHyperbolic * Simulation::method(int i) const throw()
{
  return method_list_[i];
}

//----------------------------------------------------------------------

MethodTimestep * Simulation::timestep() const throw()
{
  return timestep_;
}

//----------------------------------------------------------------------

MethodControl * Simulation::control() const throw()
{
  return control_;
}

//======================================================================

void Simulation::initialize_simulation_(Parameters * parameters) throw()
{

  //--------------------------------------------------
  // Physics parameters
  //   dimensions
  //--------------------------------------------------

  parameters->set_current_group("Physics");

  dimension_ = parameters->value_integer("dimensions",0);


  parameters->set_current_group("Domain");

  // get extent_length and check consistency

  int extent_length = parameters->list_length("extent");

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

    extent_[i] = parameters->list_value_scalar(i, "extent", extent_default);

    if (i % 2 == 1) {
      char error_message[100];
      sprintf (error_message,
	       "Parameter Domain:extent[%g] not lower than Domain:extent[%g]",
	       extent_[i-1],extent_[i]);
      ASSERT ("Simulation::Simulation", error_message, (extent_[i-1] < extent_[i]));
    }
  }
}

//----------------------------------------------------------------------

void Simulation::initialize_mesh_(Parameters * parameters) throw()
{
  parameters->set_current_group ("Mesh");

  // Block sizes

  int block_size[3];

  block_size[0] = parameters->list_value_integer(0,"block_size",4);
  block_size[1] = parameters->list_value_integer(1,"block_size",4);
  block_size[2] = parameters->list_value_integer(2,"block_size",4);

  // assume constant block size
  mesh_->set_min_block_size (block_size[0]);
  mesh_->set_max_block_size (block_size[0]);
  
  // Patch sizes

  int patch_size[3];

  patch_size[0] = parameters->list_value_integer(0,"patch_size",128);
  patch_size[1] = parameters->list_value_integer(1,"patch_size",128);
  patch_size[2] = parameters->list_value_integer(1,"patch_size",128);

  // minimum patch size is one block
  mesh_->set_min_patch_size (block_size[0]);
  mesh_->set_max_patch_size (patch_size[0]);

  // Mesh size

  int root_size[3];

  root_size[0] = parameters->list_value_integer(0,"root_size",1);
  root_size[1] = parameters->list_value_integer(1,"root_size",1);
  root_size[2] = parameters->list_value_integer(2,"root_size",1);

  mesh_->set_root_size(root_size[0],root_size[1],root_size[2]);

  parameters->set_current_group("Domain");

  double domain_lower[3];

  domain_lower[0] = parameters->list_value_scalar(0,"extent",0.0);
  domain_lower[1] = parameters->list_value_scalar(2,"extent",0.0);
  domain_lower[2] = parameters->list_value_scalar(4,"extent",0.0);

  mesh_->set_lower(domain_lower[0],
		   domain_lower[1],
		   domain_lower[2]);

  double domain_upper[3];

  domain_upper[0] = parameters->list_value_scalar(1,"extent",1.0);
  domain_upper[1] = parameters->list_value_scalar(3,"extent",1.0);
  domain_upper[2] = parameters->list_value_scalar(5,"extent",1.0);

  mesh_->set_upper(domain_upper[0],
		   domain_upper[1],
		   domain_upper[2]);

  // Parallel layout of the root mesh patch

  
  Layout * layout = new Layout;

  parameters->set_current_group("Mesh");

  int process_first = parameters->value_integer("root_process_first",0);
  int process_count = parameters->value_integer("root_process_count",1);

  layout->set_process_range(process_first, process_count);

  int block_count[3];

  block_count[0] = parameters->list_value_integer(0,"root_blocks",1);
  block_count[1] = parameters->list_value_integer(1,"root_blocks",1);
  block_count[2] = parameters->list_value_integer(2,"root_blocks",1);

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

  mesh_->set_refine    (parameters->value_integer("refine",    2));
  mesh_->set_max_level (parameters->value_integer("max_level", 0));
  mesh_->set_balanced  (parameters->value_logical("balanced",  true));
  mesh_->set_backfill  (parameters->value_logical("backfill",  true));
  mesh_->set_coalesce  (parameters->value_logical("coalesce",  true));

}

//----------------------------------------------------------------------

void Simulation::initialize_data_(Parameters * parameters) throw()
{
  FieldDescr * field_descr = data_descr_->field_descr();

  parameters->set_current_group("Field");

  // Add data fields

  int i;
  for (i=0; i<parameters->list_length("fields"); i++) {
    field_descr->insert_field
      (parameters->list_value_string(i, "fields"));
  }

  // Define default ghost zone depth for all fields, default value of 1

  int gx = parameters->list_value_integer(0,"ghosts",1);
  int gy = parameters->list_value_integer(1,"ghosts",1);
  int gz = parameters->list_value_integer(2,"ghosts",1);

  if (dimension_ < 2) gy = 0;
  if (dimension_ < 3) gz = 0;

  for (i=0; i<field_descr->field_count(); i++) {
    field_descr->set_ghosts(i,gx,gy,gz);
  }

}

//----------------------------------------------------------------------

void Simulation::initialize_control_(Parameters * parameters) throw()
{
  // Initialize stopping criteria paramaters

  parameters->set_current_group ("Stopping");

  double time_final = parameters->value_scalar("time",2.5);
  int   cycle_final = parameters->value_integer("cycle",1000);

  // Initialize output parameters

  parameters->set_current_group ("Output");

  int  cycle_dump    = parameters->value_integer("cycle_dump",10);

  // Initialize monitor parameters

  parameters->set_current_group ("Monitor");

  int  cycle_image    = parameters->value_integer("cycle_image",10);
  int  cycle_progress = parameters->value_integer("cycle_progress",1);

  // // Initial progress and image monitoring

  INCOMPLETE_MESSAGE("Simulation::initialize_control_","");

}

//----------------------------------------------------------------------

void Simulation::initialize_timestep_(Parameters * parameters) throw()
{
  INCOMPLETE_MESSAGE("Simulation::initialize_timestep_","");
}

//----------------------------------------------------------------------

void Simulation::initialize_initial_(Parameters * parameters) throw()
{
  INCOMPLETE_MESSAGE("Simulation::initialize_initial_","");
}

//----------------------------------------------------------------------

void Simulation::initialize_method_(Parameters * parameters) throw()
{
  //--------------------------------------------------
  parameters->set_current_group("Method");
  //--------------------------------------------------

  int method_count = parameters->list_length("sequence");

  if (method_count == 0) {
    ERROR_MESSAGE ("Simulation::initialize",
		   "List parameter 'Method sequence' must have length greater than zero");
  }

  for (int i=0; i<method_count; i++) {

    std::string method_name = parameters->list_value_string(i,"sequence");

    MethodHyperbolic * method = create_method_(method_name);
    if (method) method_list_.push_back(method); 

    method->initialize(data_descr_);
  }


  INCOMPLETE_MESSAGE("Simulation::initialize_method_","");
}

