// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "simulation.hpp"
#include "data.hpp" 
#include "enzo.hpp"
#include "user.hpp" 

Simulation::Simulation(Global * global)
  : domain_lower_(),
    domain_upper_(),
    global_    (global),
    mesh_      (NULL),
    user_control_(NULL),
    user_timestep_(NULL),
    user_method_(),
    data_descr_()
{

  // --------------------------------------------------
  // Initialize Domain
  // --------------------------------------------------

  // Parameter Domain::extent

  Parameters * parameters = global_->parameters();

  parameters->set_current_group("Domain");
  
  int extent_length = parameters->list_length("extent");

  if ( ! ((extent_length == 2) ||
	  (extent_length == 4) ||
	  (extent_length == 6)) ) {
    ERROR_MESSAGE ("Simulation::Simulation",
		   "'Domain extent' list must have length 2, 4, or 6");
  }

  int i;
  for (i=0; i<extent_length; i+=2) {
    domain_lower_.push_back(parameters->list_value_scalar(i,  "extent",0));
    domain_upper_.push_back(parameters->list_value_scalar(i+1,"extent",1));
    if ( ! (domain_lower_[i/2] < domain_upper_[i/2]) ) {
      ERROR_MESSAGE ("Simulation::Simulation",
		     "'domain_lower_[] not lower than domain_upper_[]");
    }
  }

  // --------------------------------------------------
  // Initialize AMR grid
  // --------------------------------------------------

  parameters->set_current_group("Mesh");

  mesh_ = new Mesh;

  // assert extent_length = 2 | 4 | 6
  mesh_->set_dimension(extent_length / 2);

  // Parameter Mesh::root_size

  std::vector<int> root_size;
  root_size.push_back(parameters->list_value_integer(0,"root_size",1));
  root_size.push_back(parameters->list_value_integer(1,"root_size",1));
  root_size.push_back(parameters->list_value_integer(2,"root_size",1));

  mesh_->set_root_size (root_size);

  // Parameter Mesh::patch_size

  mesh_->set_min_patch_size(parameters->list_value_integer(0,"patch_size",4));
  mesh_->set_max_patch_size(parameters->list_value_integer(1,"patch_size",128));

  // Parameter Mesh::block_size

  mesh_->set_min_block_size(parameters->list_value_integer(0,"block_size",4));
  mesh_->set_max_block_size(parameters->list_value_integer(1,"block_size",128));

  // Parameter Mesh::[AMR settings]

  mesh_->set_max_level     (parameters->value_integer("max_level", 0));
  mesh_->set_refine        (parameters->value_integer("refine",    2));
  mesh_->set_balanced      (parameters->value_logical("balanced",  true));
  mesh_->set_backfill      (parameters->value_logical("backfill",  true));
  mesh_->set_coalesce      (parameters->value_logical("coalesce",  true));

  // --------------------------------------------------
  // Initiazize data descriptors
  // --------------------------------------------------

  // Initialize field descriptor

  data_descr_ = new DataDescr();

  FieldDescr * field_descr = data_descr_->field_descr();

  parameters->set_current_group("Field");

  for (i=0; i<parameters->list_length("fields"); i++) {
    field_descr->insert_field
      (parameters->list_value_string(i, "fields", "unknown"));
  }

  // INITIALIZE USER METHODS

  user_initialize_();

}

//----------------------------------------------------------------------

void Simulation::initialize() throw()
{
}

//----------------------------------------------------------------------

void Simulation::advance(double stop_time, int stop_cycle) throw()
{
  INCOMPLETE_MESSAGE("Simulation::advance","Not implemented");
}
