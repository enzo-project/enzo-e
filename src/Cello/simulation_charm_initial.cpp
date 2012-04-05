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

void Block::entry_initial()
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

  // // NOTE: CHARM++ contribute() barrier is to prevent race conditions
  // // where Block::entry_refresh_face() is called before Block::entry_initial()

  // // Refresh before prepare()

  contribute( CkCallback(CkIndex_Block::entry_call_refresh(), thisProxy) );
}

//----------------------------------------------------------------------

void Simulation::entry_initial () throw()
{
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

	patch->block_array().entry_initial_enforce();

      }

    } else {

      initial->enforce(hierarchy,field_descr);

    }

  } else {

    ItPatch it_patch(hierarchy);
    Patch * patch;
    while (( patch = ++it_patch )) {
      patch->block_array().entry_call_refresh();
    }
  }
}

//----------------------------------------------------------------------

void Block::entry_initial_enforce()
{
  Simulation * simulation  = proxy_simulation.ckLocalBranch();

  Initial    * initial     = simulation->problem()->initial();
  Hierarchy  * hierarchy   = simulation->hierarchy();
  FieldDescr * field_descr = simulation->field_descr();

  initial->enforce(hierarchy,field_descr,this);
}

// //----------------------------------------------------------------------

// // Initial::init()

// //----------------------------------------------------------------------

// void Block::entry_write (int index_initial)
// {
//   Simulation * simulation = proxy_simulation.ckLocalBranch();

//   FieldDescr * field_descr = simulation->field_descr();
//   Initial * initial = simulation->problem()->initial(index_initial);

//   initial->write_block(this,field_descr,0,0,0);

//   // Synchronize via main chare before writing
//   Hierarchy * hierarchy = simulation->hierarchy();
//   int num_blocks = hierarchy->patch(0)->num_blocks();
//   simulation->proxy_block_reduce().entry_initial_reduce (num_blocks);
// }

// //----------------------------------------------------------------------

// void Block::entry_read ()
// {
//   INCOMPLETE("Block::entry_read");
// }

// //----------------------------------------------------------------------

// void BlockReduce::entry_initial_reduce(int count)
// {
//   if (++count_initial_ >= count) {
//     proxy_simulation.entry_initial_reduce();
//     count_initial_ = 0;
//   }
// }

// //----------------------------------------------------------------------

// void Simulation::entry_initial_reduce() throw()
// {
//   problem()->initial_reduce(this);
// }

// //----------------------------------------------------------------------

// void Problem::initial_reduce(Simulation * simulation) throw()
// {
//   Initial * output = initial(index_initial_);

//   output->end_write_patch();

//   int ip       = CkMyPe();
//   int ip_writer = output->process_writer();

//   if (ip != ip_writer) {

//     int n=1;  char * buffer = 0;

//     // Copy / alias buffer array of data to send
//     output->prepare_remote(&n,&buffer);

//     // Remote call to receive data
//     proxy_simulation[ip_writer].entry_initial_write (n, buffer);

//     // Close up file
//     output->close();

//     // Deallocate from prepare_remote()
//     output->cleanup_remote(&n,&buffer);

//     // Prepare for next initial
//     output->finalize();

//     // Continue with next initial object if any
//     initial_next(simulation);

//   } else {

//     proxy_simulation[ip].entry_initialial_write(0,0);

//   }

// }

// //----------------------------------------------------------------------

// void Simulation::entry_initial_write (int n, char * buffer) throw()
// {
//   problem()->initial_write(this,n,buffer);
// }

// //----------------------------------------------------------------------

// void Problem::initial_write 
// (
//  Simulation * simulation,
//  int n, char * buffer
// ) throw()
// {
//   Initial * output = initial(index_initial_);

//   if (n != 0) {
//     output->update_remote(n, buffer);
//   }

//   // fan-in from writers
//   int remaining = output->counter()->remaining();

//   if (remaining == 0) {

//     output->close();

//     output->finalize();

//     initial_next(simulation);
//   }

// }

//======================================================================

#endif /* CONFIG_USE_CHARM */


