// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 const char         parameter_file[],
 int                n,
 int                index) throw ()
  : EnzoSimulation(parameter_file, n, index)
{
#ifdef CONFIG_USE_PROJECTIONS
  traceRegisterUserEvent("Compute",10);
#endif
  // Monitor output from root pe only
  monitor_->set_active(CkMyPe() == 0);

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

  PARALLEL_PRINTF ("%s:%d DEBUG\n",__FILE__,__LINE__);

  ItPatch it_patch(hierarchy_);
  Patch * patch;
  while (( patch = ++it_patch )) {
    if (patch->blocks_allocated()) {
      patch->block_array().p_initial();
    }
  }

  PARALLEL_PRINTF ("%s:%d DEBUG\n",__FILE__,__LINE__);
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
