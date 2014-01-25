// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-01
/// @brief    Functions implementing CHARM++ output-related functions
///
/// This file contains member functions for various CHARM++ chares and
/// classes used for Output in a CHARM++ simulation.  Functions are
/// listed in roughly the order of flow-of-control.
///
///    OUTPUT
///
///    CommBlock::p_output()
///       simulation->update_state(dt,stop)
///       simulation_charm->p_output()
///
///    SimulationCharm::p_output()
///       if (block_sync_.done())
///          contribute(SimulationCharm::c_output())
///
///    SimulationCharm::c_output()
///       problem()->output_reset()   
///       problem()->output_next()   
///
///    Problem::output_next()
///       if (output)
///          output->init()
///          output->open()
///          output write_simulation()
///       else
///          simulation()->monitor_output()
///
///    Output::write_simulation() [virtual]
///       write_simulation_()
///
///    Output::write_simulation_()
///       write_hierarchy()
///
///    Output::write_hierarchy() [virtual]
///       write_hierarchy_()
///
///    Output::write_hierarchy_()
///       if (is_root())
///          block_array.p_write()
///
///    Simulation::monitor_output()
///       performance_output()
///       memory.reset_high()
///       >>>>> c_compute() >>>>>
///
///    CommBlock  p_write()
///       output.write_block(this)
///       simulation_charm.s_write()
///
///    Output::write_block() [virtual]
///       write_block_()
///
///    Output::write_block_()
///       write_field_block()
///
///    Output::write_field_block() [virtual]
///
///    Simulation::s_write()
///       if (block_sync.done())
///          contribute (SimulationCharm::c_write())
///
///    SimulationCharm::c_write()
///       problem()->output_wait()
///
///    Problem::output_wait()
///       if (ip==writer)
///          proxy_simulation[ip].p_output_write()
///       else
///          output.prepare_remote()
///          proxy_simulaion[ip_writer].p_output_write()
///          output.close()
///          output.cleanup_remote()
///          output.finalize()
///          output_next()
///
///    SimulationCharm::p_output_write()
///       problem.output_write()
///
///    Problem::output_write()
///       output().update_remote()
///       if (output->sync().done())
///          output.close()
///          output.finalize()
///          output_next()

#define TRACE_OUTPUT

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void CommBlock::p_output(CkReductionMsg * msg)
{
  
  switch_performance_ (perf_output,__FILE__,__LINE__);

  // index_.print("BEGIN PHASE OUTPUT");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  TRACE("CommBlock::p_output()");
  double * min_reduce = (double * )msg->getData();

  double dt_forest   = min_reduce[0];
  bool   stop_forest = min_reduce[1] == 1.0 ? true : false;
  set_dt   (dt_forest);
  //  printf("DEBUG CommBlock::p_output(): cycle %d dt=%f min_reduce %25.15f stop %d\n",
  //	 cycle_,dt_forest,min_reduce[1],stop_forest);

  delete msg;

  simulation->update_state(cycle_,time_,dt_forest,stop_forest);

  // Wait for all blocks to check in before calling Simulation::p_output()
  // for next output

  TRACE("CommBlock::p_output() calling SimulationCharm::p_output");
  SimulationCharm * simulation_charm = proxy_simulation.ckLocalBranch();
  simulation_charm->p_output();
}

//----------------------------------------------------------------------

void SimulationCharm::p_output ()
{
  TRACE("SimulationCharm::p_output");
  TRACE2 ("block_sync: %d/%d",block_sync_.index(),block_sync_.stop());
  if (block_sync_.next()) {
    TRACE("SimulationCharm::p_output calling c_output");
    CkCallback callback (CkIndex_SimulationCharm::c_output(), thisProxy);
    contribute(0,0,CkReduction::concat,callback);
  }
}

//----------------------------------------------------------------------

void SimulationCharm::c_output()
{
  TRACE("OUTPUT SimulationCharm::c_output()");
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
  simulation_charm->s_write();
}

//----------------------------------------------------------------------

void SimulationCharm::s_write()
{
  TRACE("SimulationCharm::s_write()");
  TRACE2 ("block_sync: %d/%d",block_sync_.index(),block_sync_.stop());
  if (block_sync_.next()) {
    CkCallback callback (CkIndex_SimulationCharm::c_write(), thisProxy);
    contribute(0,0,CkReduction::concat,callback);

  }
}

//----------------------------------------------------------------------

void SimulationCharm::c_write()
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

    proxy_simulation[ip].p_output_write(0,0);

  } else {

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

  c_compute();
}

//======================================================================


