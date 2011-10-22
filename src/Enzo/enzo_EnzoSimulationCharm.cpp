// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "simulation.hpp"
#include "enzo.hpp"

#include "simulation_charm.hpp"
#include "mesh_charm.hpp"

//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 const char         parameter_file[],
 int                n,
 CProxy_BlockReduce proxy_block_reduce, 
 int                index) throw ()
  : EnzoSimulation(parameter_file, n, proxy_block_reduce,index)
{

#ifdef CONFIG_USE_PROJECTIONS
  traceRegisterUserEvent("Compute",10);
#endif

  // Initialize on all processes
  initialize();

  // Start simulation
  run();

}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::run() throw()
{
  
  //--------------------------------------------------
  // Initial [block]
  //--------------------------------------------------

  ItPatch it_patch(hierarchy_);
  Patch * patch;
  while (( patch = ++it_patch )) {
    if (patch->blocks_allocated()) {
      patch->block_array().p_initial();
    }
  }
}


//======================================================================

#endif /* CONFIG_USE_CHARM */
