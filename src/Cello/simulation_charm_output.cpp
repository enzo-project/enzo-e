// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ output-related functions
///
/// This file contains member functions for various CHARM++ chares and
/// classes used for Output in a CHARM++ simulation.  Functions are
/// listed in roughly the order of flow-of-control.

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"

#include "simulation_charm.hpp"
#include "mesh_charm.hpp"

//----------------------------------------------------------------------

// (Called from BlockReduce::p_prepare())

void Simulation::p_output ()
{
  problem()->output_first();
  problem()->output_next(this);
}

//----------------------------------------------------------------------

void Problem::output_first() throw()
{
  index_output_ = -1;
}

//----------------------------------------------------------------------

void Problem::output_next(Simulation * simulation) throw()
{
  int cycle   = simulation->cycle();
  double time = simulation->time();

  Output * output;

  // skip over unscheduled outputs

  do {

    ++index_output_;
    output = this->output(index_output_);

  } while (output && ! output->is_scheduled(cycle, time));

  // assert ! output || output->is_scheduled(cycle_, time_)

  // output if any scheduled, else proceed with refresh

  if (output != NULL) {

    // Prepare for IO
    output->init();

    // Open files
    output->open();

    // Write hierarchy

    output->write_simulation(simulation);


  } else {

    simulation->c_monitor();

  }
}

//----------------------------------------------------------------------

// Output::init()

//----------------------------------------------------------------------

void Simulation::p_write(int index)
{
  DEBUG ("Simulation::p_write()");
  ItPatch it_patch(hierarchy_);
  Patch * patch;
  while (( patch = ++it_patch )) {
    CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
    proxy_patch->p_write(index);
  }
}

//----------------------------------------------------------------------

void Patch::p_write(int index)
{
  DEBUG ("Patch::p_write()");
  block_array_->p_write(index);
}

//----------------------------------------------------------------------

void Block::p_write (int index_output)
{
  DEBUG ("Block::p_write()");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  Output * output = simulation->problem()->output(index_output);

  output->write_block(this,field_descr,0,0,0);

  // Synchronize after writing
  proxy_patch_.s_write();
}

//----------------------------------------------------------------------

void Patch::s_write()
{
  DEBUG ("Patch::s_write()");
  if (block_counter_.remaining() == 0) {
    proxy_simulation.s_write();
  }
}

//----------------------------------------------------------------------

void Simulation::s_write()
{
  DEBUG ("Simulation::s_write()");
  problem()->output_wait(this);
}

//----------------------------------------------------------------------

void Problem::output_wait(Simulation * simulation) throw()
{
  Output * output = this->output(index_output_);

  output->end_write_patch();

  int ip       = CkMyPe();
  int ip_writer = output->process_writer();

  if (ip != ip_writer) {

    int n=1;  char * buffer = 0;

    // Copy / alias buffer array of data to send
    output->prepare_remote(&n,&buffer);

    // Remote call to receive data
    proxy_simulation[ip_writer].p_output_write (n, buffer);

    // Close up file
    output->close();

    // Deallocate from prepare_remote()
    output->cleanup_remote(&n,&buffer);

    // Prepare for next output
    output->finalize();

    // Continue with next output object if any
    output_next(simulation);

  } else {

    proxy_simulation[ip].p_output_write(0,0);

  }

}

//----------------------------------------------------------------------

void Simulation::p_output_write (int n, char * buffer)
{
  DEBUG ("Simulation::p_output_write()");
  problem()->output_write(this,n,buffer);
}

//----------------------------------------------------------------------

void Problem::output_write 
(
 Simulation * simulation,
 int n, char * buffer
) throw()
{
  DEBUG ("Problem::output_write()");
  Output * output = this->output(index_output_);

  if (n != 0) {
    output->update_remote(n, buffer);
  }

  // fan-in from writers
  int remaining = output->counter()->remaining();

  if (remaining == 0) {

    output->close();

    output->finalize();

    output_next(simulation);
  }

}

//======================================================================

#endif /* CONFIG_USE_CHARM */


