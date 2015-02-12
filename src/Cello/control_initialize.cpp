// See LICENSE_CELLO file for license and copyright information

/// @file     control_initialize.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void Simulation::initialize() throw()
{

  set_phase(phase_initial);

  initialize_config_();
  parameters_->set_monitor(false);

  initialize_monitor_();
  initialize_memory_();
  initialize_performance_();
  initialize_simulation_();

  initialize_data_descr_();

  problem_->initialize_boundary(config_,parameters_);
  problem_->initialize_initial (config_,parameters_,field_descr_);
  problem_->initialize_refine  (config_,parameters_,field_descr_);
  problem_->initialize_stopping(config_);
  problem_->initialize_output  (config_,field_descr_,factory());
  problem_->initialize_method  (config_,field_descr_);
  problem_->initialize_refresh (config_,field_descr_);
  problem_->initialize_prolong (config_);
  problem_->initialize_restrict (config_);

  initialize_hierarchy_();

  // initialize_forest_() is called in charm_initialize
  // using QD to ensure that initialize_hierarchy() is called
  // on all processors before Blocks are created

  // Barrier (why?)

  CkCallback callback 
    (CkIndex_Simulation::r_initialize_forest(NULL), thisProxy);

  // --------------------------------------------------
  contribute(0,0,CkReduction::concat,callback);
  // --------------------------------------------------
}

//----------------------------------------------------------------------

void Simulation::r_initialize_forest(CkReductionMsg * msg) 
{

  delete msg;
  
  initialize_forest_();

  // --------------------------------------------------
  // ENTRY: #2 Simulation::r_initialize_forest() -> Simulation::r_initialize_hierarchy()
  // ENTRY: callback   
  // --------------------------------------------------
  CkCallback callback 
    (CkIndex_Simulation::r_initialize_hierarchy(NULL), thisProxy);
  // --------------------------------------------------

  contribute(0,0,CkReduction::concat,callback);
}

//----------------------------------------------------------------------

void Simulation::r_initialize_hierarchy(CkReductionMsg * msg) 
{
  delete msg;

  if (CkMyPe() == 0) {
    
    // --------------------------------------------------
    // ENTRY: #3 Simulation::r_initialize_hierarchy() -> Block::p_adapt_mesh()
    // ENTRY: Block Array if Simulation is_root()
    // --------------------------------------------------
    (*hierarchy()->block_array() ).p_initial_exit();
    // --------------------------------------------------
  }
}

