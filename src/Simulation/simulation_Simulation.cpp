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
  : lower_(),
    upper_(),
    global_    (global),
    mesh_      (NULL),
    schedule_  (NULL),
    user_descr_(NULL),
    data_descr_()
{

  // --------------------------------------------------
  // Initialize Domain
  // --------------------------------------------------

  // Parameter Domain::extent

  Parameters * parameters = global_->parameters();

  parameters->set_current_group("Domain");
  
  int extent_length = parameters->list_length("extent");

  if (extent_length != 2 && 
      extent_length != 4 && 
      extent_length != 6) {
    ERROR_MESSAGE ("Simulation::initialize",
		   "List parameter 'Domain extent' must have length 2, 4, or 6");
  }

  int i;
  for (i=0; i<extent_length; i+=2) {
    lower_.push_back(parameters->list_value_scalar(i,  "extent",0));
    upper_.push_back(parameters->list_value_scalar(i+1,"extent",1));
  }

  // --------------------------------------------------
  // Initialize AMR grid
  // --------------------------------------------------

  parameters->set_current_group("Mesh");

  mesh_ = new Mesh;

  mesh_->set_dimension(extent_length / 2);

  // Parameter Mesh::root_size

  std::vector<int> root_size;

  root_size.push_back(parameters->list_value_integer(0,"root_size",1));
  root_size.push_back(parameters->list_value_integer(1,"root_size",1));
  root_size.push_back(parameters->list_value_integer(2,"root_size",1));

  mesh_->set_root_size     (root_size);
  mesh_->set_max_level     (parameters->value_integer("max_level",     0));
  mesh_->set_refine        (parameters->value_integer("refine",        2));
  mesh_->set_balanced      (parameters->value_logical("balanced",      true));
  mesh_->set_backfill      (parameters->value_logical("backfill",      true));
  mesh_->set_coalesce      (parameters->value_logical("coalesce",      true));
  mesh_->set_min_patch_size(parameters->value_integer("min_patch_size",0));
  mesh_->set_max_patch_size(parameters->value_integer("max_patch_size",0));

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

  // --------------------------------------------------
  // Initiazize user method descriptors
  // --------------------------------------------------

  parameters->set_current_group("Method");

  int method_count = parameters->list_length("sequence");

  if (method_count == 0) {
    ERROR_MESSAGE ("Simulation::initialize",
		   "List parameter 'Method sequence' must have length greater than zero");
  }

  // CREATE ENZO DESCRIPTOR OBJECT

  user_descr_ = new EnzoUserDescr(global_);

  for (i=0; i<method_count; i++) {

    std::string method_name = parameters->list_value_string(i,"sequence");

    UserMethod * user_method = user_descr_->add_user_method(method_name);

    user_method->initialize(data_descr_);
  }

  user_descr_ -> set_user_control("ignored");
  user_descr_ -> set_user_timestep("ignored");

  // --------------------------------------------------
  // Initiazize the scheduler
  // --------------------------------------------------

  schedule_ = new Schedule;
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
