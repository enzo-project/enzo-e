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

void CommBlock::r_output(CkReductionMsg * msg)
{
  
  switch_performance_ (perf_output,__FILE__,__LINE__);

  // index_.print("BEGIN PHASE OUTPUT",-1,2,false,simulation());
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  TRACE("CommBlock::r_output()");
  double * min_reduce = (double * )msg->getData();

  double dt_forest   = min_reduce[0];
  bool   stop_forest = min_reduce[1] == 1.0 ? true : false;
  set_dt   (dt_forest);
  //  printf("DEBUG CommBlock::r_output(): cycle %d dt=%f min_reduce %25.15f stop %d\n",
  //	 cycle_,dt_forest,min_reduce[1],stop_forest);

  delete msg;

  simulation->update_state(cycle_,time_,dt_forest,stop_forest);

  TRACE("CommBlock::r_output() calling SimulationCharm::output");
  SimulationCharm * simulation_charm = proxy_simulation.ckLocalBranch();

  simulation_charm->output();
}

//----------------------------------------------------------------------

void SimulationCharm::output ()
{
  TRACE("SimulationCharm::output");
  TRACE2 ("block_sync: %d/%d",block_sync_.index(),block_sync_.stop());
  if (block_sync_.next()) {
    TRACE("SimulationCharm::output calling r_output");
    CkCallback callback (CkIndex_SimulationCharm::r_output(), thisProxy);
    // --------------------------------------------------
    // ENTRY: #1 SimulationCharm::output()-> SimulationCharm::r_output()
    // ENTRY: contribute() if block_sync_.next()
    // --------------------------------------------------
    contribute(0,0,CkReduction::concat,callback);
    // --------------------------------------------------
  }
}

//----------------------------------------------------------------------

void SimulationCharm::r_output()
{
  TRACE("OUTPUT SimulationCharm::r_output()");
  problem()->output_reset();
  problem()->output_next(this);
}

//----------------------------------------------------------------------

void Problem::output_next(Simulation * simulation) throw()
{
  TRACE("OUTPUT Problem::output_next()");
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

    simulation->monitor_output();

  }
}

//----------------------------------------------------------------------

void CommBlock::p_write (int index_output)
{
  TRACE("OUTPUT CommBlock::p_write()");
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
  TRACE("SimulationCharm::write_()");
  TRACE2 ("block_sync: %d/%d",block_sync_.index(),block_sync_.stop());
  if (block_sync_.next()) {

    CkCallback callback (CkIndex_SimulationCharm::r_write(), thisProxy);

    // --------------------------------------------------
    // ENTRY: #2 SimulationCharm::write_()-> SimulationCharm::r_write()
    // ENTRY: contribute() if block_sync_.next()
    // --------------------------------------------------
    contribute(0,0,CkReduction::concat,callback);
    // --------------------------------------------------

  }
}

//----------------------------------------------------------------------

void SimulationCharm::r_write()
{
  problem()->output_wait(this);
}

//----------------------------------------------------------------------

void Problem::output_wait(Simulation * simulation) throw()
{
  TRACE("OUTPUT Problem::output_wait()");
  Output * output = this->output(index_output_);

  int ip       = CkMyPe();
  TRACE1 ("output = %p",output);
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
  TRACE("OUTPUT SimulationCharm::p_output_write()");
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
  TRACE2("OUTPUT Problem::output_write() %d %p",n,output);

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

void SimulationCharm::monitor_output()
{
  TRACE("Simulation::monitor_output()");
  Simulation::monitor_output();

  TRACE ("END   PHASE OUTPUT [SIMULATION]");

  compute();
}

//======================================================================


