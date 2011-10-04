// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ output-related functions
/// @todo      Move dt / cycle update to separate function from p_output()
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

void Simulation::p_output 
( int cycle, double time, double dt, bool stop ) throw()
{

  // Update Simulation cycle and time from reduction to main
  
  cycle_ = cycle;
  time_  = time;
  dt_    = dt;
  stop_  = stop;

  // reset output "loop" over output objects
  output_first();

  // process first output object, which continues with refresh() if done
  output_next();
}

//----------------------------------------------------------------------

void Simulation::output_first() throw()
{
  index_output_ = 0;
}

//----------------------------------------------------------------------

void Simulation::output_next() throw()
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  // find next output

  while (index_output_ < num_output() && 
	 ! output(index_output_)->is_scheduled(cycle_, time_))
    ++index_output_;

  // output if any scheduled, else proceed with refresh

  if (index_output_ < num_output()) {

    // Prepare for IO
    output(index_output_)->init();

    // Open files
    output(index_output_)->open();

    // Write hierarchy
    output(index_output_)->write_hierarchy(field_descr, hierarchy_,index_output_);


  } else {

    refresh();

  }
}

//----------------------------------------------------------------------

// Output::init()

//----------------------------------------------------------------------

void Block::p_write (int index_output)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  simulation->output(index_output)->write_block(field_descr,this,0,0,0);

  // Synchronize via main chare before writing
  int num_blocks = simulation->hierarchy()->patch(0)->num_blocks();
  simulation->proxy_block_reduce().p_output_reduce (num_blocks);
}

//----------------------------------------------------------------------

// Output::write_block()

//----------------------------------------------------------------------

void BlockReduce::p_output_reduce(int count)
{
  if (++count_output_ >= count) {
    proxy_simulation.p_output_reduce();
    count_output_ = 0;
  }
}

//----------------------------------------------------------------------

void Simulation::p_output_reduce() throw()
{
  Output * output = Simulation::output(index_output_);
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

    // Clean up and prepare for next output
    output->finalize();

    // Continue with next output object if any
    output_next();

  } else {

    proxy_simulation[ip].p_output_write(0,0);

  }

}

//----------------------------------------------------------------------

void Simulation::p_output_write (int n, char * buffer) throw()
{
  Output * output = Simulation::output(index_output_);
  int ip       = CkMyPe();
  int ip_writer = output->process_writer();

  if (n != 0) {
    output->update_remote(n, buffer);
  }

  output->counter_increment();

  int count = output->counter_value();

  if (count == output->process_stride()) {

    output->counter_reset();

    // Close up file
    output->close();

    // Clean up and prepare for next output
    output->finalize();



    output_next();
  }


}

//======================================================================

#endif /* CONFIG_USE_CHARM */


