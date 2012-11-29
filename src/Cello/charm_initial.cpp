// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_initial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-04-03
/// @brief    Functions implementing initialization functions requiring
/// CHARM++
///
/// This file contains member functions for various CHARM++ chares and
/// classes used for calling Initial objects in a CHARM++ simulation.
/// Functions are listed in roughly the order of flow-of-control.

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"


//----------------------------------------------------------------------

void SimulationCharm::initial ()
{
  problem()->initial_reset();
  problem()->initial_next(this);
}

//======================================================================

void Problem::initial_next(Simulation * simulation) throw()
{

  // find next initialization object

  Initial * initial;

  initial = this->initial(++index_initial_);

  Hierarchy *  hierarchy   = simulation->hierarchy();
  FieldDescr * field_descr = simulation->field_descr();

  if (initial != NULL) {

    if (initial->expects_blocks_allocated()) {

      DEBUG1 ("Start Initial(%d) A",index_initial_);

      ItPatch it_patch(hierarchy);
      Patch * patch;

      while (( patch = ++it_patch )) {

	CProxy_Patch * patch_proxy = (CProxy_Patch *)patch;
	patch_proxy->p_initial();

      }

    } else {

      DEBUG1 ("Start Initial(%d) B",index_initial_);

      initial->enforce_block((Block *)NULL, field_descr, hierarchy);

    }

  } else {

    SimulationCharm * simulation_charm  = 
      dynamic_cast<SimulationCharm *> (proxy_simulation.ckLocalBranch());

    simulation_charm->refresh();

  }
}

//----------------------------------------------------------------------

void Patch::p_initial()
{
  TRACE("Patch::p_initial()");
  TRACE("Patch::p_initial(Patch) calling Block::p_initial()");
  block_array()->p_initial();
}

//----------------------------------------------------------------------

void Patch::s_initial()
{
  if (block_loop_.done()) {
    proxy_simulation.s_initial();
  }
}

//----------------------------------------------------------------------

void SimulationCharm::s_initial()
{
  if (patch_loop_.done()) {
    delete parameters_;
    parameters_ = 0;
    p_refresh();
  }
}

//======================================================================

#endif /* CONFIG_USE_CHARM */


