// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodCheckpoint.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-05-14
/// @brief    Implementation of the Checkpoint "method"

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

#define TRACE_METHOD_CHECKPOINT(BLOCK,MSG) CkPrintf ("TRACE_METHOD_CHECKPOINT %d %s: %s\n", \
                                                     CkMyPe(),BLOCK->name().c_str(),MSG); fflush(stdout);

MethodCheckpoint::MethodCheckpoint
( std::vector< std::string > path_name)
  : Method(),
    path_name_(path_name)
{
  CkPrintf ("MethodCheckpoint::MethodCheckpoint\n");
  Refresh * refresh = cello::refresh(ir_post_);

  refresh->add_all_fields();
}
//----------------------------------------------------------------------

void MethodCheckpoint::pup (PUP::er &p)
{
  TRACEPUP;
  Method::pup(p);
  p | path_name_;
}

//----------------------------------------------------------------------

void MethodCheckpoint::compute ( Block * block) throw()
{
  // barrier to ensure tree traversal doesn't reach a block
  // before the method has started
  TRACE_METHOD_CHECKPOINT(block,"compute");
  CkCallback callback(CkIndex_Block::r_method_checkpoint_continue(nullptr), 
                      cello::block_array());
  block->contribute(callback);

}

//----------------------------------------------------------------------

void Block::r_method_checkpoint_continue(CkReductionMsg *msg)
{
  TRACE_METHOD_CHECKPOINT(this,"r_method_checkpoint_compute");
  delete msg;
  MethodCheckpoint * method = static_cast<MethodCheckpoint*> (this->method());
  method->compute_continue(this);
}

//----------------------------------------------------------------------

void MethodCheckpoint::compute_continue(Block * block)
{
  TRACE_METHOD_CHECKPOINT(block,"compute_continue");
  if (block->index().is_root()) {
    static int counter = 0;
    std::string path_name = cello::directory(&path_name_,counter++,block);
    proxy_main.p_checkpoint_method(CkNumPes(),path_name);
  } else {
    block->compute_done();
  }
}

//----------------------------------------------------------------------

void Simulation::r_write_checkpoint_method()
{
  // ERROR need to return to MethodCheckpoint with block containing is_root() somehow
  // save block proxy maybe
  ERROR("Simulation::r_write_checkpoint_method()",
        "Not implemented yet");
  performance_->start_region(perf_output);
  create_checkpoint_link();
  problem()->output_wait(this);
  performance_->stop_region(perf_output);
}

