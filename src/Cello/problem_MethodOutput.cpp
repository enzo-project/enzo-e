// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOutput.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-27
/// @brief    Implementation of the Output "method"

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

#define TRACE_METHOD_OUTPUT(block,msg)                                                \
  {                                                                     \
    CkPrintf ("TRACE_METHOD_OUTPUT %s %s\n",msg,block->name().c_str()); \
    fflush(stdout);                                                     \
  }


//----------------------------------------------------------------------

MethodOutput::MethodOutput
(std::vector<int> field_list,
 std::vector<int> particle_list,
 int ghost_depth,
 int min_face_rank,
 bool all_fields,
 bool all_particles,
 int blocking_x,
 int blocking_y,
 int blocking_z)
  : Method(),
    field_list_(field_list),
    particle_list_(particle_list),
    ghost_depth_(ghost_depth),
    min_face_rank_(min_face_rank),
    all_fields_(all_fields),
    all_particles_(all_particles),
    blocking_()
{
  Refresh * refresh = cello::refresh(ir_post_);

  if (ghost_depth > 0) refresh->set_ghost_depth (ghost_depth);
  refresh->set_min_face_rank(min_face_rank);

  refresh->add_field("density");
  // add fields to refresh
  if (all_fields_) {
    refresh->add_all_fields();
  } else {
    const int nf=field_list_.size();
    for (int i_f=0; i_f<nf; i_f++) {
      cello::refresh(ir_post_)->add_field(field_list_[i_f]);
    }
  }

  // add particles to refresh
  if (all_particles_) {
    refresh->add_all_particles();
  } else {
    const int nf=particle_list_.size();
    for (int i_f=0; i_f<nf; i_f++) {
      cello::refresh(ir_post_)->add_particle(particle_list_[i_f]);
    }
  }

  // initialize blocking partition
  blocking_[0] = blocking_x;
  blocking_[1] = blocking_y;
  blocking_[2] = blocking_z;
}

//----------------------------------------------------------------------

void MethodOutput::compute ( Block * block) throw()
{
  TRACE_METHOD_OUTPUT(block,"compute");

  // barrier to ensure tree traversal doesn't reach a block
  // before the method has started

  CkCallback callback(CkIndex_Block::r_method_output_continue(nullptr), 
                      cello::block_array());
  block->contribute(callback);

}

//----------------------------------------------------------------------

void Block::r_method_output_continue(CkReductionMsg *msg)
{
  delete msg;
  MethodOutput * method = static_cast<MethodOutput*> (this->method());
  method->compute_continue(this);
}

//----------------------------------------------------------------------

void MethodOutput::compute_continue(Block * block)
{
  TRACE_METHOD_OUTPUT(block,"compute_continue");

  if (is_writer_(block->index())) {
    int a3[3];
    block->index().array(a3,a3+1,a3+2);
    int index_min[3] = {a3[0] - (a3[0] % blocking_[0]),
                        a3[1] - (a3[1] % blocking_[1]),
                        a3[2] - (a3[2] % blocking_[2])};
    int index_max[3] = {index_min[0] + blocking_[0],
                        index_min[1] + blocking_[1],
                        index_min[2] + blocking_[2]};
    BlockTrace bt (cello::rank(),index_min,index_max);
    MsgOutput * msg_output = new MsgOutput(bt,this);
    msg_output->set_index_send (block->index());

    CkPrintf ("%s: File Open\n",block->name().c_str());

    if (block->is_leaf()) {

      CkPrintf ("%s: File Write %s\n",
                block->name().c_str(),
                block->name().c_str());
    }

    BlockTrace * block_trace = msg_output->block_trace();
    
    if (!block_trace->next(block->is_leaf())) {
      Index index_next = block_trace->top();
      CkPrintf ("%s: Send to %s\n",
                block->name().c_str(),
                block->name(index_next).c_str());
      cello::block_array()[index_next].p_method_output_next(msg_output);
      
    } else {
      delete msg_output;
      CkPrintf ("%s: File Close\n",block->name().c_str());
      TRACE_METHOD_OUTPUT(block,"done");
      block->compute_done();
    }
  }

  if (block->level() < 0) {
    // negative level blocks do not participate
    TRACE_METHOD_OUTPUT(block,"done");
    block->compute_done();
  }
}

//----------------------------------------------------------------------

void Block::p_method_output_next (MsgOutput * msg_output)
{
  msg_output->method_output_->next(this,msg_output);
}

//----------------------------------------------------------------------

void MethodOutput::next(Block * block, MsgOutput * msg_output )
{
  TRACE_METHOD_OUTPUT(block,"next");
  const bool is_leaf = block->is_leaf();
  BlockTrace * bt = msg_output->block_trace();
  bt->next(block->is_leaf());
  msg_output->set_index_send (block->index());
  if (is_leaf) {
    // if leaf, copy data to msg and send to writer
    Index index_home = bt->home();
    cello::block_array()[index_home].p_method_output_write(msg_output);
  } else {
    // if non-leaf, forward to child
    Index index_next = bt->top();
    cello::block_array()[index_next].p_method_output_next(msg_output);
  }
  TRACE_METHOD_OUTPUT(block,"done");
  block->compute_done();

}

//----------------------------------------------------------------------

void Block::p_method_output_write (MsgOutput * msg_output)
{
  msg_output->method_output_->write(this,msg_output);
}

//----------------------------------------------------------------------

void MethodOutput::write(Block * block, MsgOutput * msg_output )
{
  TRACE_METHOD_OUTPUT(block,"write");
  CkPrintf ("%s: File Write %s\n",
            block->name().c_str(),
            block->name(msg_output->index_send()).c_str());
  // delete data
  msg_output->set_index_send (block->index());
  BlockTrace * bt = msg_output->block_trace();
  Index index_next = bt->top();
  Index index_home = bt->home();
  if (index_next == index_home) {
    // done
    CkPrintf ("%s: File Close\n",block->name().c_str());
    block->compute_done();
  } else {
    cello::block_array()[index_next].p_method_output_next(msg_output);
  }  
}

//----------------------------------------------------------------------

int MethodOutput::is_writer_ (Index index) 
{
  int a3[3];
  index.array(a3,a3+1,a3+2);
  int level = index.level();
  return ((level==0) &&
          ( a3[0] % blocking_[0] == 0) &&
          ( a3[1] % blocking_[1] == 0) &&
          ( a3[2] % blocking_[2] == 0));
}

