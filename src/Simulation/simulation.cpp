// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      
/// @brief     
 
//====================================================================
// PUBLIC FUNCTIONS
//====================================================================

Simulation::Simulation()
  : lower_(),
    upper_(),
    amr_(NULL),
    parallel_(NULL),
    fields_(NULL),
    methods_(NULL)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a Simulation object
 *
 *********************************************************************
 */
{
}

Simulation::~Simulation()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Delete a Simulation object
 *
 *********************************************************************
 */
{
}
    
/// Initialize from parameters
void Simulation::initialize (Parameters * parameters)
{
  // --------------------------------------------------
  // Initialize domain
  // --------------------------------------------------

  parameters->set_group("Domain");
  
  int extent_length = parameters->list_length("extent");
  if (extent_length % 2 != 0) {
    ERROR_MESSAGE ("Simulation::initialize",
		   "List parameter 'Domain extent' should have an even length");
  }
  assert (extent_length % 2 == 0);
  for (int i=0; i<extent_length/2; i++) {
    lower_.push_back(parameters->list_value_integer(i));
    upper_.push_back(parameters->list_value_integer(i+1));
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
  // Initialize parallel object
  // --------------------------------------------------

  parallel_ = new Parallel ();

  // --------------------------------------------------
  // Initialize data fields
  // --------------------------------------------------

  fields_(NULL),

  // --------------------------------------------------
      // 
  // --------------------------------------------------
    methods_(NULL)
}

//====================================================================
// PRIVATE FUNCTIONS
//====================================================================

