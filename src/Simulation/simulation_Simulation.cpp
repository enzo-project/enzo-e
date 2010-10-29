// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "simulation.hpp"
#include "data.hpp" 
#include "user.hpp" 

Simulation::Simulation(Global * global)
  : dimension_(0),
    domain_lower_(),
    domain_upper_(),
    global_    (global),
    mesh_      (NULL),
    user_control_(NULL),
    user_timestep_(NULL),
    user_method_(),
    data_descr_()
{

  Parameters * parameters = global_->parameters();

  // Initialize dimension_

  parameters->set_current_group("Physics");

  dimension_ = parameters->value_integer("dimensions",0);

  // Initialize domain extents

  parameters->set_current_group("Domain");

  // get extent_length and check consistency

  int extent_length = parameters->list_length("extent");

  bool valid_extent_length = ((extent_length == 2) ||
			      (extent_length == 4) ||
			      (extent_length == 6));
  ASSERT ("Simulation::Simulation",
	  "Parameter Domain:extent list must have length 2, 4, or 6",
	  valid_extent_length);

  ASSERT ("Simulation::Simulation",
	  "Parameter mismatch between Physics:dimensions and Domain:extent",
	  (dimension_ == extent_length / 2));

  // initialize domain_lower_[] and domain_upper_[]

  for (int i=0; i<extent_length; i+=2) {

    int k = i/2;

    domain_lower_[k] = parameters->list_value_scalar(i,  "extent",0.0);
    domain_upper_[k] = parameters->list_value_scalar(i+1,"extent",1.0);

    ASSERT ("Simulation::Simulation",
	    "Parameter Domain:lower not lower than Domain:upper",
	    (domain_lower_[k] < domain_upper_[k]));
  }

//   // assert extent_length = 2 | 4 | 6
//   mesh_->set_dimension(extent_length / 2);

  // Parameter Mesh::root_size

//   std::vector<int> root_size;
//   root_size.push_back(parameters->list_value_integer(0,"root_size",1));
//   root_size.push_back(parameters->list_value_integer(1,"root_size",1));
//   root_size.push_back(parameters->list_value_integer(2,"root_size",1));

//   mesh_->set_root_size (root_size);

  // Parameter Mesh::patch_size

  // --------------------------------------------------
  // Initiazize data descriptors
  // --------------------------------------------------

  // Initialize field descriptor

  data_descr_ = new DataDescr;

  FieldDescr * field_descr = data_descr_->field_descr();

  parameters->set_current_group("Field");

  for (int i=0; i<parameters->list_length("fields"); i++) {
    field_descr->insert_field
      (parameters->list_value_string(i, "fields", "unknown"));
  }

  // --------------------------------------------------
  // Initialize AMR grid
  // --------------------------------------------------

  initialize_mesh_();
}

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
{
  delete mesh_;
  delete user_control_;
  delete user_timestep_;
  delete data_descr_;
}

//----------------------------------------------------------------------
Simulation::Simulation(const Simulation & classname) throw()
{
  INCOMPLETE_MESSAGE("Simulation::Simulation","");
}

//----------------------------------------------------------------------

Simulation & Simulation::operator= (const Simulation & classname) throw()
{
  INCOMPLETE_MESSAGE("Simulation::operator =","");
  return *this;
}

//----------------------------------------------------------------------

int Simulation::dimension() const throw()
{
  return dimension_; 
};

//----------------------------------------------------------------------

int Simulation::domain (int domain_lower[], int domain_upper[]) throw()
{
  if (dimension_ >= 1) {
    domain_lower[0] = domain_lower_[0];
    domain_upper[0] = domain_upper_[0];
  }
  if (dimension_ >= 2) {
    domain_lower[1] = domain_lower_[1];
    domain_upper[1] = domain_upper_[1];
  }
  if (dimension_ >= 3) {
    domain_lower[2] = domain_lower_[2];
    domain_upper[2] = domain_upper_[2];
  }
  return dimension_;
}

//----------------------------------------------------------------------

void Simulation::initialize_mesh_() throw()
{

  mesh_ = new Mesh(data_descr_);

  Parameters * parameters = global_->parameters();

  parameters->set_current_group ("Mesh");

  // Patch size range

  int min_patch_size = parameters->list_value_integer(0,"patch_size",4);
  int max_patch_size = parameters->list_value_integer(1,"patch_size",128);

  mesh_->set_min_patch_size (min_patch_size);
  mesh_->set_max_patch_size (max_patch_size);

  // Block size range

  int min_block_size = parameters->list_value_integer(0,"block_size",4);
  int max_block_size = parameters->list_value_integer(1,"block_size",128);
  
  mesh_->set_min_block_size (min_block_size);
  mesh_->set_max_block_size (max_block_size);
  
  // Mesh size

  int root_size[3];

  root_size[0] = parameters->list_value_integer(0,"root_size",1);
  root_size[1] = parameters->list_value_integer(1,"root_size",1);
  root_size[2] = parameters->list_value_integer(2,"root_size",1);

  mesh_->set_root_size(root_size[0],root_size[1],root_size[2]);

  // Lower and upper domain extents

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

    // Create the root patch

  mesh_->set_root_patch(new Patch (data_descr_,
				   root_size,
				   max_block_size,
				   domain_lower,domain_upper));

  // Mesh AMR settings

  mesh_->set_refine    (parameters->value_integer("refine",    2));
  mesh_->set_max_level (parameters->value_integer("max_level", 0));
  mesh_->set_balanced  (parameters->value_logical("balanced",  true));
  mesh_->set_backfill  (parameters->value_logical("backfill",  true));
  mesh_->set_coalesce  (parameters->value_logical("coalesce",  true));

}
//----------------------------------------------------------------------

void Simulation::advance(double stop_time, int stop_cycle) throw()
{
  INCOMPLETE_MESSAGE("Simulation::advance","Not implemented");
}
