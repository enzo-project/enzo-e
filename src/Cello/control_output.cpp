// See LICENSE_CELLO file for license and copyright information

/// @file     control_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ output-related functions
/// @ingroup  Control

/* #define TRACE_OUTPUT */

#ifdef TRACE_OUTPUT
#  define TRACE_LOCAL(M) printf ("TRACE Output %s:%d" M "\n", __FILE__,__LINE__); fflush(stdout);
#else
#  define TRACE_LOCAL(M) /*  */
#endif

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void CommBlock::output_begin_ ()
{

  simulation()->set_phase(phase_output);

  TRACE_LOCAL("CommBlock::output_begin_()");

  int cycle   = simulation()->cycle();
  double time = simulation()->time();

  Output * output;

  // skip over unscheduled outputs

  int index_output = -1;

  do {
    output = simulation()->problem()->output(++index_output);

  } while (output && ! output->is_scheduled(cycle, time));

  // invariant: (output == NULL) || output->is_scheduled()

  TRACE2("output_begin %d %p",index_output,output);

  if (output != NULL) {

    ((SimulationCharm * )simulation()) -> begin_output();

  } else {

    output_exit_();

  }
}

//----------------------------------------------------------------------

void SimulationCharm::begin_output ()
{

  TRACE_LOCAL("CommBlock::output_begin()");

#ifdef CONFIG_USE_MEMORY
  trace_mem_ = Memory::instance()->bytes() - trace_mem_;
#endif

  performance()->switch_region(perf_output,__FILE__,__LINE__);

  if (sync_output_begin_.next()) {

    CkCallback callback (CkIndex_SimulationCharm::r_output(NULL), thisProxy);
    // --------------------------------------------------
    // ENTRY: #1 SimulationCharm::output()-> SimulationCharm::r_output()
    // ENTRY: contribute() if sync_output_begin_.next()
    // --------------------------------------------------
    contribute(0,0,CkReduction::concat,callback);
    // --------------------------------------------------
  }
}

//----------------------------------------------------------------------

void SimulationCharm::r_output(CkReductionMsg * msg)
{
 
  TRACE_LOCAL("SimulationCharm::r_output()");

  delete msg;

  problem()->output_reset();
  problem()->output_next(this);
}

//----------------------------------------------------------------------

void Problem::output_next(Simulation * simulation) throw()
{

  simulation->set_phase(phase_output);

  TRACE_LOCAL("Problem::output_next()");

  int cycle   = simulation->cycle();
  double time = simulation->time();

  Output * output;

  // skip over unscheduled outputs

  do {

    output = this->output(++index_output_);

  } while (output && ! output->is_scheduled(cycle, time));

  // output if any scheduled, else proceed with refresh

  if (output != NULL) {

    output->init();
    output->open();
    output->write_simulation(simulation);
    output->next();

  } else {

    ((SimulationCharm *)simulation)->output_exit();

  }
}

//----------------------------------------------------------------------

void CommBlock::p_output_write (int index_output)
{
  TRACE_LOCAL("CommBlock::p_output_write()");
  FieldDescr * field_descr = simulation()->field_descr();
  Output * output = simulation()->problem()->output(index_output);

  output->write_block(this,field_descr);

  ((SimulationCharm *) simulation())->write_();
}

//----------------------------------------------------------------------

void SimulationCharm::write_()
{
  TRACE_LOCAL("SimulationCharm::write_()");
  if (sync_output_write_.next()) {

    CkCallback callback (CkIndex_SimulationCharm::r_write(NULL), thisProxy);

    // --------------------------------------------------
    // ENTRY: #2 SimulationCharm::write_()-> SimulationCharm::r_write()
    // ENTRY: contribute() if sync_output_write_.next()
    // --------------------------------------------------
    contribute(0,0,CkReduction::concat,callback);
    // --------------------------------------------------

  }
}

//----------------------------------------------------------------------

void SimulationCharm::r_write(CkReductionMsg * msg)
{
  TRACE_LOCAL("SimulationCharm::r_write()");
  delete msg;

  problem()->output_wait(this);
}

//----------------------------------------------------------------------

void SimulationCharm::r_write_checkpoint()
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
    // ENTRY: #3 Problem::output_wait()-> SimulationCharm::p_output_write()
    // ENTRY: writer Simulation if is writer
    // --------------------------------------------------
    proxy_simulation[ip].p_output_write(0,0);
    // --------------------------------------------------

  } else {

    int n=0;  char * buffer = 0;

    // Copy / alias buffer array of data to send
    output->prepare_remote(&n,&buffer);

    // Remote call to receive data

    // --------------------------------------------------
    // ENTRY: #4 Problem::output_wait()-> SimulationCharm::p_output_write()
    // ENTRY: writer Simulation if not writer
    // --------------------------------------------------
    proxy_simulation[ip_writer].p_output_write (n, buffer);
    // --------------------------------------------------

    // Close up file

    output->close();

    // Deallocate from prepare_remote()
    output->cleanup_remote(&n,&buffer);

    // Prepare for next output
    output->finalize();

    // Continue with next output object if any
    output_next(simulation);

  }

}

//----------------------------------------------------------------------

void SimulationCharm::p_output_write (int n, char * buffer)
{
  TRACE_LOCAL("SimulationCharm::p_output_write()");
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

  if (output->sync_write()->next()) {

    output->close();

    output->finalize();

    output_next(simulation);
  }

}

//----------------------------------------------------------------------

void SimulationCharm::output_exit()
{
  // reset debug output files to limit file size
  debug_close();
  debug_open();

  TRACE_LOCAL("SimulationCharm::output_exit()");
  if (hierarchy()->group_process()->is_root()) 
    hierarchy()->block_array()->p_output_exit();
}

//======================================================================


