// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ output-related functions

#define TRACE_OUTPUT

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void SimulationCharm::begin_output ()
{
#ifdef TRACE_MEMORY
  trace_mem_ = Memory::instance()->bytes() - trace_mem_;
#endif

  performance()->switch_region(perf_output,__FILE__,__LINE__);

  if (block_sync_.next()) {

    CkCallback callback (CkIndex_SimulationCharm::r_output(NULL), thisProxy);
    // --------------------------------------------------
    // ENTRY: #1 SimulationCharm::output()-> SimulationCharm::r_output()
    // ENTRY: contribute() if block_sync_.next()
    // --------------------------------------------------
    contribute(0,0,CkReduction::concat,callback);
    // --------------------------------------------------
  }
}

//----------------------------------------------------------------------

void SimulationCharm::r_output(CkReductionMsg * msg)
{

  delete msg;

  problem()->output_reset();
  problem()->output_next(this);
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

  // output if any scheduled, else proceed with refresh

  if (output != NULL) {

    output->init();
    output->open();
    output->write_simulation(simulation);

  } else {

    ((SimulationCharm *)simulation)->output_exit();

  }
}

//----------------------------------------------------------------------

void CommBlock::p_output_write (int index_output)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  Output * output = simulation->problem()->output(index_output);

  output->write_block(this,field_descr);

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();

  simulation_charm->write_();
}

//----------------------------------------------------------------------

void SimulationCharm::write_()
{
  if (block_sync_.next()) {

    CkCallback callback (CkIndex_SimulationCharm::r_write(NULL), thisProxy);

    // --------------------------------------------------
    // ENTRY: #2 SimulationCharm::write_()-> SimulationCharm::r_write()
    // ENTRY: contribute() if block_sync_.next()
    // --------------------------------------------------
    contribute(0,0,CkReduction::concat,callback);
    // --------------------------------------------------

  }
}

//----------------------------------------------------------------------

void SimulationCharm::r_write(CkReductionMsg * msg)
{
  delete msg;

  problem()->output_wait(this);
}

//----------------------------------------------------------------------

void Problem::output_wait(Simulation * simulation) throw()
{
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

    int n=1;  char * buffer = 0;

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
  problem()->output_write(this,n,buffer); 
}

//----------------------------------------------------------------------

void Problem::output_write 
(
 Simulation * simulation,
 int n, char * buffer
) throw()
{
  Output * output = this->output(index_output_);

  if (n != 0) {
    output->update_remote(n, buffer);
  }

  if (output->sync()->next()) {

    output->close();

    output->finalize();

    output_next(simulation);
  }

}

//----------------------------------------------------------------------

void SimulationCharm::output_exit()
{
  TRACE("Simulation::monitor_output()");

  monitor_-> print("", "-------------------------------------");
  monitor_-> print("Simulation", "cycle %04d", cycle_);
  monitor_-> print("Simulation", "time-sim %15.12f",time_);
  monitor_-> print("Simulation", "dt %15.12g", dt_);

  performance_output ();

  Memory::instance()->reset_high();

#ifdef TRACE_MEMORY
  trace_mem_ = Memory::instance()->bytes() - trace_mem_;
  PARALLEL_PRINTF ("memory output %lld\n",trace_mem_);
#endif

#ifdef TRACE_MEMORY
  trace_mem_ = Memory::instance()->bytes();
#endif

#ifdef CELLO_DEBUG
  fprintf (fp_debug(),"%s:%d cycle %03d\n", __FILE__,__LINE__,cycle_);
#endif

  if (cycle_ > 0 ) {
    performance()->stop_region (perf_cycle,__FILE__,__LINE__);
  }

  if (stop_) {
    
    performance_write();

    // --------------------------------------------------
    // ENTRY: #1 SimulationCharm::compute()-> Main::p_exit()
    // ENTRY: Main if stop
    // --------------------------------------------------
    proxy_main.p_exit(CkNumPes());
    // --------------------------------------------------

  } else {

    performance()->start_region (perf_cycle,__FILE__,__LINE__);
    performance()->switch_region (perf_compute,__FILE__,__LINE__);

    if (hierarchy()->group_process()->is_root()) 
      
      // --------------------------------------------------
      // ENTRY: #2 SimulationCharm::compute()-> CommBlock::p_compute_enter()
      // ENTRY: Block Array if Simulation is_root()
      // --------------------------------------------------
      hierarchy()->block_array()->p_compute_enter(cycle_,time_,dt_);
      // --------------------------------------------------
  }
#ifdef TRACE_MEMORY
  trace_mem_ = Memory::instance()->bytes() - trace_mem_;
  PARALLEL_PRINTF ("memory compute %lld\n",trace_mem_);
#endif
}


//======================================================================


