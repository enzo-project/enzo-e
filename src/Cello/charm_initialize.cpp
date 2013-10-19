// See LICENSE_CELLO file for license and copyright information

/// @file     charm_initialize.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization

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

  Simulation::initialize();

  CkCallback callback (CkIndex_SimulationCharm::r_initialize_forest(), thisProxy);
  contribute(0,0,CkReduction::concat,callback);

}

//----------------------------------------------------------------------

void SimulationCharm::r_initialize_forest() 
{

  initialize_forest_();

  CkCallback callback (CkIndex_SimulationCharm::r_initialize_end(), thisProxy);
  contribute(0,0,CkReduction::concat,callback);
}

//----------------------------------------------------------------------

void SimulationCharm::r_initialize_end() 
{
  if (group_process_->is_root()) {
#ifdef TEMP_NEW_ADAPT
    (*hierarchy()->block_array() ).p_adapt_mesh();
#else /* TEMP_NEW_ADAPT */
    (*hierarchy()->block_array() ).p_adapt_begin();
#endif /* TEMP_NEW_ADAPT */
  }
}

//----------------------------------------------------------------------
