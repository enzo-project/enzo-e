// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "enzo.hpp"

#include "simulation_charm.hpp"
#include "mesh_charm.hpp"

#include "simulation.hpp"


//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 const char         parameter_file[],
 int                n) throw ()
  : EnzoSimulation(parameter_file, n)
{

#ifdef CONFIG_USE_PROJECTIONS
  traceRegisterUserEvent("Compute",10);
#endif

  initialize();

}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::run() throw()
{
  DEBUG("EnzoSimulationCharm::run()");
  ItPatch it_patch(hierarchy_);
  Patch * patch;

  // count patches for Patch::p_initial()
  int patch_count = 0;

  DEBUG("Counting patches");
  while (( patch = ++it_patch )) {
    // count local patches
    ++patch_count;

  }
  DEBUG1("Patch count = %d",patch_count);
    
  // set patch counter for s_patch() synchronization
  patch_counter_.set_max(patch_count + 1);

  // Initialize hierarchy

  DEBUG("Calling Patch::p_initial() loop");
  while (( patch = ++it_patch )) {
    CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
    DEBUG1("Calling %p Patch::p_initial()",proxy_patch);
    proxy_patch->p_initial();
  }
  DEBUG0;
  s_initial();
  DEBUG0;
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
