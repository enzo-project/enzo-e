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
#include "mesh.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"


//----------------------------------------------------------------------

void SimulationCharm::c_initial()
{
  ItPatch it_patch(hierarchy_);
  Patch * patch;

  // count patches for Patch::p_initial()
  int patch_count = 0;

  while (( patch = ++it_patch )) {
    // count local patches
    ++patch_count;

  }

  // set patch counter for s_patch() synchronization
  patch_counter_.set_max(patch_count + 1);

  DEBUG1 ("patch count = %d",patch_count);

  // Initialize hierarchy

  while (( patch = ++it_patch )) {
    CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
    proxy_patch->p_initial();
  }

  s_initial();

}

//----------------------------------------------------------------------

void SimulationCharm::x_request_patch
(
 int patch_id,
 int block_rank
)
{
  // call x_send_patch()
  CkChareID patch_proxy;
  proxy_simulation[block_rank].x_send_patch(patch_id,patch_proxy);
  DEBUG ("x_request_patch()");
}

//----------------------------------------------------------------------

void SimulationCharm::x_send_patch
(
 int patch_id,
 CkChareID proxy_patch
 )
{
  DEBUG ("x_send_patch()");
}

//----------------------------------------------------------------------

void SimulationCharm::p_initial ()
{
  DEBUG("SimulationCharm::p_initial()");
  // reset initial "loop" over initial objects
  problem()->initial_first();

  // process first initial object, which continues with refresh() if done
  problem()->initial_next(this);
}

//----------------------------------------------------------------------

void Problem::initial_first() throw()
{
  DEBUG("Problem::initial_first()");
  index_initial_ = 0;
}

//----------------------------------------------------------------------

void Problem::initial_next(Simulation * simulation) throw()
{
  DEBUG("Problem::initial_next()");
  // find next initial

  Initial * initial = this->initial(index_initial_);

  // initial if any scheduled, else proceed with refresh

  Hierarchy * hierarchy = simulation->hierarchy();
  FieldDescr * field_descr = simulation->field_descr();

  if (initial != NULL) {

    DEBUG1 ("Initial expects blocks allocated = %s",
	    initial->expects_blocks_allocated() ? "true" : "false");
    if (initial->expects_blocks_allocated()) {

      ItPatch it_patch(hierarchy);
      Patch * patch;

      while (( patch = ++it_patch )) {

	CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
	DEBUG("Problem::initial_next() calling Patch::p_initial()");
	proxy_patch->p_initial();

      }

    } else {

      DEBUG("Problem::initial_next() calling Initial::enforce(Hierarchy)");
      initial->enforce(hierarchy,field_descr);

    }

  } else {

    ItPatch it_patch(hierarchy);
    Patch * patch;
    while (( patch = ++it_patch )) {
      CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
      DEBUG("Problem::initial_next() calling Patch::p_refresh()");
      proxy_patch->p_refresh();
    }
  }
}

//----------------------------------------------------------------------

void Patch::p_initial()
{
  DEBUG("Patch::p_initial()");
  DEBUG1("block_array() = %p",block_array());
  block_array()->p_initial();
}

//----------------------------------------------------------------------

void Block::p_initial()
{
  DEBUG("Block::p_initial()");
  Simulation * simulation  = proxy_simulation.ckLocalBranch();
  DEBUG0;
  FieldDescr * field_descr = simulation->field_descr();

  // Initialize the block

  DEBUG0;
  allocate(field_descr);
  DEBUG0;

  // Set the Block cycle and time to match Simulation's

  set_cycle(simulation->cycle());
  set_time (simulation->time());
  set_dt   (simulation->dt());

  DEBUG0;
  // Perform any additional initialization for derived class 

  initialize ();
  DEBUG0;

  // Apply the initial conditions 

  Initial * initial = simulation->problem()->initial();

  initial->enforce(simulation->hierarchy(),field_descr,this);
  DEBUG0;

  // Continue with Patch::s_initial

  proxy_patch_.s_initial();
  DEBUG0;

}

//----------------------------------------------------------------------

void Patch::s_initial()
{
  if (block_counter_.remaining() == 0) {
    DEBUG("Patch::s_initial() calling SimulationCharm::s_initial()");
    proxy_simulation.s_initial();
  } else  DEBUG("Patch::s_initial() skipping");

}

//----------------------------------------------------------------------

void Block::p_read ()
{
  INCOMPLETE("Block::p_read");
}

//======================================================================

#endif /* CONFIG_USE_CHARM */


