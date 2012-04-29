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

#ifdef NEW_INITIAL

void SimulationCharm::c_initial ()
{
  problem()->initial_reset();
  problem()->initial_next(this);
}

#endif

//======================================================================

#ifndef NEW_INITIAL

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

  // set patch loop counter for s_patch() synchronization
  patch_loop_.set_max(patch_count + 1);

  // Initialize hierarchy

  while (( patch = ++it_patch )) {
    CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
    proxy_patch->p_initial();
  }

  s_initial();

}

#endif

//----------------------------------------------------------------------

void Problem::initial_next(Simulation * simulation) throw()
{
#ifdef NEW_INITIAL

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

	CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
	proxy_patch->p_initial();

      }

    } else {

    DEBUG1 ("Start Initial(%d) B",index_initial_);

      initial->enforce(hierarchy,field_descr);

    }

  } else {

    DEBUG1 ("Start Initial(%d) C",index_initial_);

    ItPatch it_patch(hierarchy);
    Patch * patch;
    while (( patch = ++it_patch )) {
      CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
      proxy_patch->p_refresh();
    }
  }
#endif
}
//----------------------------------------------------------------------

void SimulationCharm::x_request_patch
(
 int patch_id,
 int block_rank
)
{
#ifdef NEW_INITIAL
  // call x_send_patch()
  CkChareID patch_proxy;
  proxy_simulation[block_rank].x_send_patch(patch_id,patch_proxy);
  DEBUG ("x_request_patch()");
#endif
}

//----------------------------------------------------------------------

void SimulationCharm::x_send_patch
(
 int patch_id,
 CkChareID proxy_patch
 )
{
#ifdef NEW_INITIAL

#endif
}

//----------------------------------------------------------------------

void Patch::p_initial()
{
  block_array()->p_initial();
}

//----------------------------------------------------------------------

void Block::p_initial()
{
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

  // Continue with Patch::s_initial

  proxy_patch_.s_initial();

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
    c_refresh();
  }
}

//======================================================================

void Block::p_read (int index_initial)
{
  INCOMPLETE("Block::p_read");
}
#endif /* CONFIG_USE_CHARM */


