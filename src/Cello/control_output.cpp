// See LICENSE_CELLO file for license and copyright information

/// @file     control_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @ingroup  Control
/// @brief    Functions implementing CHARM++ output-related functions

// #define TRACE_OUTPUT

#ifdef TRACE_OUTPUT
#  define TRACE_LOCAL(M) printf ("TRACE Output %s:%d %d" M "\n", __FILE__,__LINE__,CkMyPe()); fflush(stdout);
#else
#  define TRACE_LOCAL(M) /*  */
#endif

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void Block::output_begin_ ()
{
  TRACE_LOCAL("Block::output_begin_()");

  // Determine if there is any output this cycle

  simulation()->set_phase(phase_output);

  int cycle   = simulation()->cycle();
  double time = simulation()->time();

  Output * output;

  int index_output = -1;

  do {

    output = simulation()->problem()->output(++index_output);

  } while (output && ! output->is_scheduled(cycle, time));

  if (output != NULL) {

    // Start output if any...

    simulation() -> begin_output();

  } else {

    // ...otherwise continue with next phase

    output_exit_();

  }
}

//----------------------------------------------------------------------

void Simulation::begin_output ()
{

  TRACE_LOCAL("Block::output_begin()");

  // Switching from Block to Simulation: wait for last Block

  if (sync_output_begin_.next()) {

    performance()->switch_region(perf_output,__FILE__,__LINE__);

    // Barrier

    TRACE_LOCAL("Block::output_begin() calling Simulation::r_output()");
    // --------------------------------------------------
    //    CkCallback callback (CkIndex_Simulation::r_output(NULL), thisProxy);
    //    contribute(0,0,CkReduction::concat,callback);
    Simulation * simulation = proxy_simulation.ckLocalBranch();
    simulation->r_output(NULL);
    // --------------------------------------------------
  }
}

//----------------------------------------------------------------------

void Simulation::r_output(CkReductionMsg * msg)
{
 
  TRACE_LOCAL("Simulation::r_output()");

  delete msg;

  // start first output

  problem()->output_reset();
  problem()->output_next(this);
}

//----------------------------------------------------------------------

void Problem::output_next(Simulation * simulation) throw()
{

  TRACE_LOCAL("Problem::output_next()");

  simulation->set_phase(phase_output);

  int cycle   = simulation->cycle();
  double time = simulation->time();

  Output * output;

  // Find next schedule output

  do {

    output = this->output(++index_output_);

  } while (output && ! output->is_scheduled(cycle, time));

  if (output != NULL) {

    // Perform output if any...

    output->init();
    output->open();
    output->write_simulation(simulation);
    output->next();

  } else {

    // ...otherwise exit output phase

    simulation->output_exit();

  }
}

//----------------------------------------------------------------------

void Block::p_output_write (int index_output)
{
  TRACE_LOCAL("Block::p_output_write()");

  FieldDescr * field_descr = simulation()->field_descr();
  Output * output = simulation()->problem()->output(index_output);

  output->write_block(this,field_descr);

  simulation()->write_();
}

//----------------------------------------------------------------------

void Simulation::write_()
{
  TRACE_LOCAL("Simulation::write_()");
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
  TRACE_LOCAL("Simulation::r_write()");
  delete msg;

  problem()->output_wait(this);
}

//----------------------------------------------------------------------

void Simulation::r_write_checkpoint()
{
  problem()->output_wait(this);
}

//----------------------------------------------------------------------

void Problem::output_wait(Simulation * simulation) throw()
{
  TRACE_LOCAL("Problem::output_wait()");
  
  Output * output = this->output(index_output_);

  int ip       = CkMyPe();
  int ip_writer = output->process_writer();

  if (ip == ip_writer) {

    // --------------------------------------------------
    proxy_simulation[ip].p_output_write(0,0);
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

  }

}

//----------------------------------------------------------------------

void Simulation::p_output_write (int n, char * buffer)
{
  TRACE_LOCAL("Simulation::p_output_write()");
  problem()->output_write(this,n,buffer); 
}

//----------------------------------------------------------------------

void Problem::output_write 
(
 Simulation * simulation,
 int n, char * buffer
) throw()
{
  TRACE_LOCAL("Problem::output_write()");

  Output * output = this->output(index_output_);

  if (n != 0) {
    output->update_remote(n, buffer);
  }

  // ERROR HERE ON RESTART WITH DIFFERENT +p
  if (output->sync_write()->next()) {
    output->close();
    output->finalize();
    output_next(simulation);
  }

}

//----------------------------------------------------------------------

void Simulation::output_exit()
{
  TRACE_LOCAL("Simulation::output_exit()");

  // reset debug output files to limit file size
  debug_close();
  debug_open();

  if (CkMyPe() == 0) hierarchy()->block_array()->p_output_end();
}

//----------------------------------------------------------------------

void Block::p_output_end()
{
  TRACE_LOCAL("Simulation_output_end()");
  control_sync(CkIndex_Block::r_stopping_enter(NULL),"contribute");
}
//======================================================================


