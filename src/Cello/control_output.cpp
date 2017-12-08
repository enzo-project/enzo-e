// See LICENSE_CELLO file for license and copyright information

/// @file     control_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @ingroup  Control
/// @brief    Functions implementing CHARM++ output-related functions

// #define TRACE_OUTPUT

#ifdef TRACE_OUTPUT
#  define TRACE_OUTPUT(M) printf ("%d TRACE_OUTPUT %s:%d " M "\n",CkMyPe(), __FILE__,__LINE__); fflush(stdout);
#else
#  define TRACE_OUTPUT(M) /*  */
#endif

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void Block::output_begin_ ()
{
  simulation() -> begin_output();
}

//----------------------------------------------------------------------

void Simulation::begin_output ()
{

  TRACE_OUTPUT("Simulation::begin_output()");

  // Switching from Block to Simulation: wait for last Block

  if (sync_output_begin_.next()) {

    set_phase(phase_output);

    performance_->start_region(perf_output);

    problem()->output_reset();
    problem()->output_next(this);

    performance_->stop_region(perf_output);
  }
}

//----------------------------------------------------------------------

void Problem::output_next(Simulation * simulation) throw()
{

  TRACE_OUTPUT("Problem::output_next()");

  simulation->set_phase(phase_output);

  int cycle   = simulation->cycle();
  double time = simulation->time();

  Output * output;

  // Find next schedule output (index_output_ initialized to -1)

  do {

    output = this->output(++index_output_);

  } while (output && ! output->is_scheduled(cycle, time));

  // assert (! output) || ( output->is_scheduled() )
  
  if (output != NULL) {

    // Perform output if any...

    const int stride = output->stride_wait();

    if (CkMyPe() % stride == 0) {

      simulation->do_output (index_output_);

    }

  } else {

    // ...otherwise exit output phase

    simulation->output_exit();

  }
}

//----------------------------------------------------------------------

void Simulation::do_output(int index_output)
{
  TRACE_OUTPUT("Simulation::p_do_output()");
  Output * output = problem()->output(index_output);
  output->init();
  output->open();
  output->write_simulation(this);
  output->next();
}

//----------------------------------------------------------------------

void Block::p_output_write (int index_output)
{
  performance_start_ (perf_output);
  
  TRACE_OUTPUT("Block::p_output_write()");

  FieldDescr * field_descr = simulation()->field_descr();
  ParticleDescr * particle_descr = simulation()->particle_descr();
  Output * output = simulation()->problem()->output(index_output);

  output->write_block(this,field_descr,particle_descr);

  simulation()->write_();
  performance_stop_ (perf_output);
}

//----------------------------------------------------------------------

void Simulation::write_()
{
  TRACE_OUTPUT("Simulation::write_()");
  if (sync_output_write_.next()) {

    // --------------------------------------------------
    //    CkCallback callback (CkIndex_Simulation::r_write(NULL), thisProxy);
    //    contribute(0,0,CkReduction::concat,callback);
    r_write(NULL);
    // --------------------------------------------------

  }
}

//----------------------------------------------------------------------

void Simulation::r_write(CkReductionMsg * msg)
{
  performance_->start_region(perf_output);
  TRACE_OUTPUT("Simulation::r_write()");
  delete msg;

  problem()->output_wait(this);
  performance_->stop_region(perf_output);
}

//----------------------------------------------------------------------

void Simulation::r_write_checkpoint()
{
  performance_->start_region(perf_output);
  TRACE_OUTPUT("Simulation::r_write_checkpoint()");
  create_checkpoint_link();
  problem()->output_wait(this);
  performance_->stop_region(perf_output);
}

//----------------------------------------------------------------------

void Problem::output_wait(Simulation * simulation) throw()
{
  TRACE_OUTPUT("Problem::output_wait()");
  
  Output * output = this->output(index_output_);

  int ip        = CkMyPe();
  int ip_writer = output->process_writer();

  if (ip == ip_writer) {

    // --------------------------------------------------
    //    proxy_simulation[ip].p_output_write(0,0);
    output_write(simulation,0,0);
    // --------------------------------------------------

  } else {

    int n=0;  char * buffer = 0;

    // Copy / alias buffer array of data to send
    output->prepare_remote(&n,&buffer);

    // Remote call to receive data
    // --------------------------------------------------
    proxy_simulation[ip_writer].p_output_write (n, buffer);
    // --------------------------------------------------

    output->close();
    output->cleanup_remote(&n,&buffer);
    output->finalize();
    output_next(simulation);
    const int stride = output->stride_wait();
    if ((CkMyPe()+1)%stride != 0 && (CkMyPe()+1) < CkNumPes()) {
      proxy_simulation[CkMyPe()+1].p_do_output(index_output_);
    }

  }

}

//----------------------------------------------------------------------

void Simulation::p_output_write (int n, char * buffer)
{
  TRACE_OUTPUT("Simulation::p_output_write()");
  problem()->output_write(this,n,buffer); 
}

//----------------------------------------------------------------------

void Problem::output_write 
(
 Simulation * simulation,
 int n, char * buffer
) throw()
{
  TRACE_OUTPUT("Problem::output_write()");

  Output * output = this->output(index_output_);

  if (n != 0) {
    output->update_remote(n, buffer);
  }

  // ERROR HERE ON RESTART WITH DIFFERENT +p
  if (output->sync_write()->next()) {
    output->close();
    output->finalize();
    output_next(simulation);
    const int stride = output->stride_wait();
    if ((CkMyPe()+1)%stride != 0 && (CkMyPe()+1) < CkNumPes()) {
      proxy_simulation[CkMyPe()+1].p_do_output(index_output_);
    }
  }

}

//----------------------------------------------------------------------

void Simulation::output_exit()
{
  TRACE_OUTPUT("Simulation::output_exit()");

  // reset debug output files to limit file size
  debug_close();
  debug_open();

  if (CkMyPe() == 0) hierarchy()->block_array()->p_output_end();
}

//----------------------------------------------------------------------

void Block::p_output_end()
{
  performance_start_(perf_output);
  TRACE_OUTPUT("Block::p_output_end()");
  performance_stop_(perf_output);
  output_exit_();
  // control_sync(CkIndex_Block::r_stopping_enter(NULL),sync_barrier);
}
//======================================================================


