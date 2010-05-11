// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "error.hpp" 
#include "monitor.hpp"
#include "simulation.hpp"
#include "parameters.hpp"

Simulation::Simulation()
  : lower_(),
    upper_(),
    mesh_(NULL),
    data_(),
    methods_(NULL)
{
}

/// Initialize from parameters
void Simulation::initialize ()
{

  Parameters * parameters = Parameters::instance();

  // --------------------------------------------------
  // Initialize Domain
  // --------------------------------------------------

  // Parameter Domain::extent

  parameters->set_current_group("Domain");
  
  int extent_length = parameters->list_length("extent");

  if (extent_length != 2 && 
      extent_length != 4 && 
      extent_length != 6) {
    ERROR_MESSAGE ("Simulation::initialize",
		   "List parameter 'Domain extent' must have length 2, 4, or 6");
  }

  for (int i=0; i<extent_length; i+=2) {
    lower_.push_back(parameters->list_value_scalar(i,  "extent",0));
    upper_.push_back(parameters->list_value_scalar(i+1,"extent",1));
  }

  // --------------------------------------------------
  // Initialize AMR grid
  // --------------------------------------------------

  parameters->set_current_group("Mesh");

  mesh_ = new Mesh();

  mesh_->set_dimension(extent_length / 2);

  // Parameter Mesh::root_size

  std::vector<int> root_size;

  root_size.push_back(parameters->list_value_integer(0,"root_size",1));
  root_size.push_back(parameters->list_value_integer(1,"root_size",1));
  root_size.push_back(parameters->list_value_integer(2,"root_size",1));

  mesh_->set_root_size(root_size);

  // Mesh::max_level

  mesh_->set_max_level(parameters->value_integer("max_level",0));

  // Mesh::refine

  mesh_->set_refine(parameters->value_integer("refine",2));

  // Mesh::balanced

  mesh_->set_balanced(parameters->value_logical("balanced",true));

  // Mesh::backfill

  mesh_->set_backfill(parameters->value_logical("backfill",true));

  // Mesh::coalesce

  mesh_->set_coalesce(parameters->value_logical("coalesce",true));

  // Mesh::min_patch_size

  mesh_->set_min_patch_size(parameters->value_integer("min_patch_size",0));

  // Mesh::max_patch_size

  mesh_->set_max_patch_size(parameters->value_integer("max_patch_size",0));

  // --------------------------------------------------
  // Initiazize methods
  // --------------------------------------------------

  parameters->set_current_group("Method");

  int num_methods = parameters->subgroup_count();
  for (int i=0; i<num_methods; i++) {
    

  }
  INCOMPLETE_MESSAGE("Simulation::initialize","Initializing Methods");

  // --------------------------------------------------
  // Initialize Mesh
  // --------------------------------------------------

  INCOMPLETE_MESSAGE("Simulation::initialize","Initializing Mesh");
  
  // --------------------------------------------------
  // Initialize data 
  // --------------------------------------------------

  INCOMPLETE_MESSAGE("Simulation::initialize","Initializing Fields");

}
