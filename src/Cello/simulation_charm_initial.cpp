// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_initial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-04-03
/// @brief    Functions implementing CHARM++ initialization-related functions
///
/// This file contains member functions for various CHARM++ chares and
/// classes used for calling Initial objects in a CHARM++ simulation.
/// Functions are listed in roughly the order of flow-of-control.

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "simulation_charm.hpp"
#include "mesh.hpp"
#include "mesh_charm.hpp"

//----------------------------------------------------------------------

void Block::p_initial()
{
  TRACE0;
  Simulation * simulation  = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();

  // Initialize the block

  allocate(field_descr);

  // Set the Block cycle and time to match Simulation's

  set_cycle(simulation->cycle());
  set_time (simulation->time());
  set_dt   (simulation->dt());

  // Perform any additional initialization for derived class 

  initialize ();

  // Apply the initial conditions 

  Initial * initial = simulation->problem()->initial();

  initial->enforce(simulation->hierarchy(),field_descr,this);

  // // NOTE: CHARM++ contribute() barrier is to prevent race conditions
  // // where Block::p_refresh_face() is called before Block::p_initial()

  // // Refresh before prepare()

  contribute( CkCallback(CkIndex_Block::p_call_refresh(), thisProxy) );
}

//----------------------------------------------------------------------

void Simulation::p_initial () throw()
{
  TRACE0;
  // reset initial "loop" over initial objects
  problem()->initial_first();

  // process first initial object, which continues with refresh() if done
  problem()->initial_next(this);
}

//----------------------------------------------------------------------

void Problem::initial_first() throw()
{
  index_initial_ = 0;
}

//----------------------------------------------------------------------

void Problem::initial_next(Simulation * simulation) throw()
{
  // find next initial

  Initial * initial = this->initial(index_initial_);

  // initial if any scheduled, else proceed with refresh

  Hierarchy * hierarchy = simulation->hierarchy();
  FieldDescr * field_descr = simulation->field_descr();

  if (initial != NULL) {

    if (initial->expects_blocks_allocated()) {

      ItPatch it_patch(hierarchy);
      Patch * patch;

      while (( patch = ++it_patch )) {

	patch->block_array().p_initial_enforce();

      }

    } else {

      initial->enforce(hierarchy,field_descr);

    }

  } else {

    ItPatch it_patch(hierarchy);
    Patch * patch;
    while (( patch = ++it_patch )) {
      patch->block_array().p_call_refresh();
    }
  }
}

//----------------------------------------------------------------------

void Block::p_initial_enforce()
{
  Simulation * simulation  = proxy_simulation.ckLocalBranch();

  Initial    * initial     = simulation->problem()->initial();
  Hierarchy  * hierarchy   = simulation->hierarchy();
  FieldDescr * field_descr = simulation->field_descr();

  initial->enforce(hierarchy,field_descr,this);
}

//======================================================================

#endif /* CONFIG_USE_CHARM */


