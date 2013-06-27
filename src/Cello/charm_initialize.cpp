// See LICENSE_CELLO file for license and copyright information

/// @file     charm_initialize.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void SimulationCharm::p_initialize_begin() 
{
  initialize();  // virtual: calls EnzoSimulationCharm::initialize()
}

//----------------------------------------------------------------------

void SimulationCharm::initialize() throw()
{
  TRACE("SimulationCharm::initialize()");

  Simulation::initialize();

  if (group_process_->is_root()) {
    CkStartQD 
      ( CkCallback(CkIndex_SimulationCharm::q_initialize_forest(), thisProxy));
  }

}

//----------------------------------------------------------------------

void SimulationCharm::q_initialize_forest() 
{
  initialize_forest_();

  if (group_process_->is_root()) {
    CkStartQD 
      ( CkCallback(CkIndex_SimulationCharm::q_initialize_end(), thisProxy));
  
  }
}

//----------------------------------------------------------------------

void SimulationCharm::q_initialize_end() 
{
  if (group_process_->is_root()) {
    (*hierarchy()->block_array() ).p_adapt_begin();
  }
}

//----------------------------------------------------------------------

#endif /* CONFIG_USE_CHARM */
