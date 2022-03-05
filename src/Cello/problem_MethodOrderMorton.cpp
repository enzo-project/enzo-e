// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOrderMorton.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief

#include "problem.hpp"

//----------------------------------------------------------------------

MethodOrderMorton::MethodOrderMorton(int min_level) throw ()
  : Method(),
    is_index_(-1),
    is_weight_(-1),
    is_weight_child_(-1),
    min_level_(min_level)
{
  Refresh * refresh = cello::refresh(ir_post_);
  cello::simulation()->refresh_set_name(ir_post_,name());
  refresh->add_field("density");

  /// Create Scalar data for ordering index
  is_index_       = cello::scalar_descr_int()->new_value(name() + ":index");
  is_count_   = cello::scalar_descr_int()->new_value(name() + ":count");
  is_weight_      = cello::scalar_descr_int()->new_value(name() + ":weight");
  is_weight_child_ = cello::scalar_descr_int()->new_value(name() + ":weight_child",
                                                          cello::num_children());
  is_sync_index_  = cello::scalar_descr_sync()->new_value(name() + ":sync_index");
  is_sync_weight_ = cello::scalar_descr_sync()->new_value(name() + ":sync_weight");
}

//======================================================================

void MethodOrderMorton::compute (Block * block) throw()
{
  // Initialize counters, then barrier to ensure counters initialized
  // before first entry method can arrive

  *pindex_(block) = 0;
  *pcount_(block) = 0;
  *pweight_(block) = 1;
  for (int i=0; i<cello::num_children(); i++) {
    *pweight_child_(block,i) = 0;
  }
  psync_index_(block)->set_stop(1 + 1);
  psync_weight_(block)->set_stop(1 + cello::num_children());
  CkCallback callback (CkIndex_Block::r_method_order_morton_continue(nullptr),
                       block->proxy_array());
  block->contribute (callback);

}

//----------------------------------------------------------------------

void Block::r_method_order_morton_continue(CkReductionMsg * msg)
{
  delete msg;
  static_cast<MethodOrderMorton*>
    (this->method())->compute_continue(this);
}

//----------------------------------------------------------------------

void MethodOrderMorton::compute_continue(Block * block)
{
  send_weight(block, 0, true);
}

//======================================================================
void MethodOrderMorton::send_weight(Block * block, int weight_child, bool self)
{
  // update own weight
  // if not at finest level, send weight to parent
  int weight = *pweight_(block);
  int ic3[3] = {0,0,0};
  if (self) {
    recv_weight(block,ic3,0,true);
  }
  const int level = block->level();
  if ((!self || block->is_leaf()) && level > min_level_)  {
    const Index index_parent = block->index().index_parent(min_level_);
    block->index().child(level,ic3,ic3+1,ic3+2,min_level_);
    cello::block_array()[index_parent].p_method_order_morton_weight(ic3,weight,block->index());
    send_index(block, 0, 0, self);
  } else if (level == min_level_) {
    send_index(block, 0, weight, self);
    if (!self) {
      CkCallback callback (CkIndex_Block::r_method_order_morton_complete(nullptr),
                           block->proxy_array());
      block->contribute (callback);
    }
  }
}

//----------------------------------------------------------------------

void Block::p_method_order_morton_weight(int ic3[3], int weight, Index index_child)
{
  static_cast<MethodOrderMorton*>
    (this->method())->recv_weight(this, ic3,weight,false);
}

//----------------------------------------------------------------------

void MethodOrderMorton::recv_weight
(Block * block, int ic3[3], int weight, bool self)
{
  // Update children weight if needed
  if (!self) {
    *pweight_(block) += weight;
    int i = ic3[0] + 2*(ic3[1]+2*ic3[2]);
    *pweight_child_(block,i) = weight;
  }
  if (psync_weight_(block)->next()) {
    // Forward weight to parent when computed
    int ic3[3] = {0,0,0};
    block->index().child(block->level(),ic3,ic3+1,ic3+2,min_level_);
    send_weight(block,*pweight_(block),false);
  }
}

void MethodOrderMorton::send_index
(Block * block, int index_parent, int count, bool self)
{
  *pcount_(block) = count;
  if (!block->is_leaf()) {
    int index = *pindex_(block) + 1;
    for (int ic=0; ic<cello::num_children(); ic++) {
      int ic3[3];
      ic3[0] = (ic>>0) & 1;
      ic3[1] = (ic>>1) & 1;
      ic3[2] = (ic>>2) & 1;
      Index index_child = block->index().index_child(ic3,min_level_);
      cello::block_array()[index_child].p_method_order_morton_index(index,count);
      index += *pweight_child_(block,ic);
    }
  }
}

void Block::p_method_order_morton_index(int index, int count)
{
  static_cast<MethodOrderMorton*>
    (this->method())->recv_index(this, index, count, false);
}

void MethodOrderMorton::recv_index
(Block * block, int index, int count, bool self)
{
  if (!self) {
    *pindex_(block)     = index;
    *pcount_(block) = count;
  }
  if (psync_index_(block)->next()) {
    send_index(block,index, count, false);
    CkCallback callback (CkIndex_Block::r_method_order_morton_complete(nullptr),
                       block->proxy_array());
    block->contribute (callback);
  }
}

//----------------------------------------------------------------------

void Block::r_method_order_morton_complete(CkReductionMsg * msg)
{
  delete msg;
  static_cast<MethodOrderMorton*>
    (this->method())->compute_complete(this);
}

//----------------------------------------------------------------------

void MethodOrderMorton::compute_complete(Block * block)
{
  block->compute_done();
}


//======================================================================

int * MethodOrderMorton::pindex_(Block * block)
{
  Scalar<int> scalar(cello::scalar_descr_int(),
                     block->data()->scalar_data_int());
  return scalar.value(is_index_);
}

//----------------------------------------------------------------------

int * MethodOrderMorton::pcount_(Block * block)
{
  Scalar<int> scalar(cello::scalar_descr_int(),
                     block->data()->scalar_data_int());
  return scalar.value(is_count_);
}

//----------------------------------------------------------------------

Sync * MethodOrderMorton::psync_index_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return scalar.value(is_sync_index_);
}

//----------------------------------------------------------------------

int * MethodOrderMorton::pweight_(Block * block)
{
  Scalar<int> scalar(cello::scalar_descr_int(),
                     block->data()->scalar_data_int());
  return scalar.value(is_weight_);
}

//----------------------------------------------------------------------

int * MethodOrderMorton::pweight_child_(Block * block, int i)
{
  Scalar<int> scalar(cello::scalar_descr_int(),
                     block->data()->scalar_data_int());
  return scalar.value(is_weight_child_)+i;
}

//----------------------------------------------------------------------

Sync * MethodOrderMorton::psync_weight_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return scalar.value(is_sync_weight_);
}

