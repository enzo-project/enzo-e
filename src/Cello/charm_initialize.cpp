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

void SimulationCharm::p_initialize() 
{
  initialize();  // virtual: calls EnzoSimulationCharm::initialize()
}

//----------------------------------------------------------------------

void SimulationCharm::initialize() throw()
{
  TRACE("SimulationCharm::initialize()");
  Simulation::initialize();

  s_initialize();
}

//----------------------------------------------------------------------

void SimulationCharm::s_initialize()
{
  TRACE("SimulationCharm::s_initialize()");

  if (group_process_->is_root()) {
    CkStartQD (CkCallback(CkIndex_CommBlock::p_phase_adapt(),*hierarchy()->block_array()));
   }

}

#endif /* CONFIG_USE_CHARM */
