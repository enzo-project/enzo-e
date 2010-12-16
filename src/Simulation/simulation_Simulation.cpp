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
    data_(0),
    control_(0),
    timestep_(0),
    initialize_list_(),
    method_list_()
{
  // Initialize parameter defaults

  dimension_ = 1;

  extent_[0] = 0.0;
  extent_[1] = 1.0;
  extent_[2] = 0.0;
  extent_[3] = 1.0;
  extent_[4] = 0.0;
  extent_[5] = 1.0;

  // Initialize simulation component defaults

  data_ = new DataDescr;

  mesh_ = new Mesh(data_);

}

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
{
  delete data_;
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

  initialize_parameters_(parameters);

  // Initialize simulation components

  initialize_mesh_(parameters);
  initialize_data_(parameters);
  initialize_control_(parameters);
  initialize_timestep_(parameters);
  initialize_initial_(parameters);
  initialize_methods_(parameters);

}

//----------------------------------------------------------------------

void Simulation::execute() throw()
{
  INCOMPLETE_MESSAGE("Simulation::execute","");
}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  INCOMPLETE_MESSAGE("Simulation::finalize","");
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

void Simulation::initialize_parameters_(Parameters * parameters) throw()
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
  ASSERT ("Simulation::initialize_parameters_",
	  "Parameter Domain:extent list must have length 2, 4, or 6",
	  valid_extent_length);

  ASSERT ("Simulation::initialize_parameters_",
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

  root_patch->set_data_descr(data_);
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
  FieldDescr * field_descr = data_->field_descr();

  parameters->set_current_group("Field");

  for (int i=0; i<parameters->list_length("fields"); i++) {
    field_descr->insert_field
      (parameters->list_value_string(i, "fields", "unknown"));
  }
}

//----------------------------------------------------------------------

void Simulation::initialize_control_(Parameters * parameters) throw()
{
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

void Simulation::initialize_methods_(Parameters * parameters) throw()
{
  INCOMPLETE_MESSAGE("Simulation::initialize_methods_","");
}

//----------------------------------------------------------------------

MethodHyperbolic * Simulation::add_method (std::string method_name) throw()
{
  MethodHyperbolic * method = create_method_(method_name);
  if (method) method_list_.push_back(method); 
  return method;
};
