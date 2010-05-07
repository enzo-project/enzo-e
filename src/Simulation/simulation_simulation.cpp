// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      
/// @brief     
 
#include "simulation.hpp"

Simulation::Simulation()
  : lower_(),
    upper_(),
    amr_(NULL),
    data_(),
    methods_(NULL)
{
}

/// Initialize from parameters
void Simulation::initialize ()
{

  Parameters * parameters = Parameters::instance();

  // --------------------------------------------------
  // Initialize domain
  // --------------------------------------------------

  parameters->set_group("Domain");
  
  int extent_length = parameters->list_length("extent");
  if (extent_length==2 || extent_length==4 || extent_length==6) {
    ERROR_MESSAGE ("Simulation::initialize",
		   "List parameter 'Domain extent' must have length 2, 4, or 6");
  }

  for (int i=0; i<extent_length/2; i++) {
    lower_.push_back(parameters->list_value_integer(i,  "extent",0));
    upper_.push_back(parameters->list_value_integer(i+1,"extent",1));
  }

  // --------------------------------------------------
  // Initialize AMR grid
  // --------------------------------------------------

  parameters->set_group("Grid");

  amr_ = new Amr();

  
  //  root      = [400,400]; # size of the root grid
  //  levels    = 4;         # maximum effective levels (for r_factor = 2)
  //  refine    = 4;         # refinement factor
  //  type      = tree;      # AMR type: patch or tree
  //  full_tree = false;     # whether tree is full (ala Flash) or not
  //  backfill  = true;      # whether to backfill for refinement > 2
  //  coalesce  = true;      # whether to coalesce small patches to one big one
  //  patch_min = 4;         # minimum patch size
  //  patch_max = 128;       # maximum patch size

  // --------------------------------------------------
  // Initiazize methods
  // --------------------------------------------------

  INCOMPLETE_MESSAGE("Simulation::initialize","Initializing Methods");

  // --------------------------------------------------
  // Initialize Amr / Arrays
  // --------------------------------------------------

  INCOMPLETE_MESSAGE("Simulation::initialize","Initializing Amr");
  
  // --------------------------------------------------
  // Initialize data 
  // --------------------------------------------------

  INCOMPLETE_MESSAGE("Simulation::initialize","Initializing Fields");

}
