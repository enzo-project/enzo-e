// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
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
 CProxy_BlockReduce proxy_block_reduce) throw ()
  : EnzoSimulation(parameter_file, n, proxy_block_reduce)
{

#ifdef CONFIG_USE_PROJECTIONS
  traceRegisterUserEvent("Compute",10);
#endif

  initialize();

  run();

}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::run() throw()
{
  
  // Call Block::entry_initial() on all blocks

  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while (( patch = ++it_patch )) {

    if (patch->blocks_allocated()) {

      patch->block_array().entry_initial();

    } else {

      // Initial      

      // Blocks don't exist: read Blocks from file and insert into Patch

      // UNCOMMENTING THE FOLLOWING EXITS PREMATURELY FOR NON-RESTART
       // Initial * initial = problem()->initial();
       // initial->enforce(hierarchy(),field_descr());
    }

  }
}


//======================================================================

#endif /* CONFIG_USE_CHARM */
