// See LICENSE_CELLO file for license and copyright information

/// @file     control_initialize.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void SimulationCharm::initialize() throw()
{

  set_phase(phase_initial);

  Simulation::initialize();

  // --------------------------------------------------
  // ENTRY: #1 SimulationCharm::initialize() -> SimulationCharm::r_initialize_forest()
  // ENTRY: callback
  // --------------------------------------------------
  CkCallback callback 
    (CkIndex_SimulationCharm::r_initialize_forest(NULL), thisProxy);
  contribute(0,0,CkReduction::concat,callback);
  // --------------------------------------------------
}

//----------------------------------------------------------------------

void SimulationCharm::r_initialize_forest(CkReductionMsg * msg) 
{

  delete msg;
  
  initialize_forest_();

  // --------------------------------------------------
  // ENTRY: #2 SimulationCharm::r_initialize_forest() -> SimulationCharm::r_initialize_hierarchy()
  // ENTRY: callback   
  // --------------------------------------------------
  CkCallback callback 
    (CkIndex_SimulationCharm::r_initialize_hierarchy(NULL), thisProxy);
  // --------------------------------------------------

  contribute(0,0,CkReduction::concat,callback);
}

//----------------------------------------------------------------------

void SimulationCharm::r_initialize_hierarchy(CkReductionMsg * msg) 
{
  delete msg;

  if (group_process_->is_root()) {
    
    // --------------------------------------------------
    // ENTRY: #3 SimulationCharm::r_initialize_hierarchy() -> CommBlock::p_adapt_mesh()
    // ENTRY: Block Array if Simulation is_root()
    // --------------------------------------------------
    (*hierarchy()->block_array() ).p_adapt_enter();
    // --------------------------------------------------
  }
}

//----------------------------------------------------------------------
