// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOrderMorton.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-01
/// @brief

#include "problem.hpp"

// #define TRACE_ORDER

#ifdef TRACE_ORDER
#  define TRACE_ORDER_BLOCK(MSG,BLOCK)          \
  CkPrintf ("TRACE_ORDER %s %s\n",              \
            std::string(MSG).c_str(),           \
            BLOCK->name().c_str());             \
  fflush(stdout);
#else
#  define TRACE_ORDER_BLOCK(MSG,BLOCK) /* ... */
#endif


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
  const int n = cello::num_children();
  is_index_        = cello::scalar_descr_long_long()->new_value(name() + ":index");
  is_count_        = cello::scalar_descr_long_long()->new_value(name() + ":count");
  is_next_         = cello::scalar_descr_index()->new_value(name() + ":next");
  is_weight_       = cello::scalar_descr_long_long()->new_value(name() + ":weight");
  is_weight_child_ = cello::scalar_descr_long_long()->new_value(name() + ":weight_child",n);
  is_sync_index_  = cello::scalar_descr_sync()->new_value(name() + ":sync_index");
  is_sync_weight_ = cello::scalar_descr_sync()->new_value(name() + ":sync_weight");
}

//======================================================================

void MethodOrderMorton::compute (Block * block) throw()
{
  // Initialize counters, then barrier to ensure counters initialized
  // before first entry method can arrive
  TRACE_ORDER_BLOCK("compute",block);
  Sync * sync_index = psync_index_(block);
  Sync * sync_weight = psync_weight_(block);

  *pindex_(block) = 0;
  *pcount_(block) = 0;
  *pweight_(block) = 1;
  for (int i=0; i<cello::num_children(); i++) {
    *pweight_child_(block,i) = 0;
  }
  sync_index->reset();
  sync_weight->reset();
  sync_index->set_stop(1 + 1);
  sync_weight->set_stop(1 + cello::num_children());

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
  TRACE_ORDER_BLOCK("continue",block);
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
    TRACE_ORDER_BLOCK("send_weight",block);
    cello::block_array()[index_parent].p_method_order_morton_weight
      (ic3,weight,block->index());
    send_index(block, 0, 0, self);
  } else if (level == min_level_) {

    const int rank = cello::rank();
    int na3[3];
    cello::simulation()->hierarchy()->root_blocks(na3,na3+1,na3+2);
    Index index_next = block->index().next
      (rank,na3,block->is_leaf(),min_level_);

    *pindex_(block) = 0;
    *pcount_(block) = 0;
    *pnext_(block) = index_next;

    send_index(block, 0, weight, self);
    if (!self) {
      CkCallback callback
        (CkIndex_Block::r_method_order_morton_complete (nullptr),
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
  TRACE_ORDER_BLOCK("recv_weight",block);
  // Update children weight if needed
  if (!self) {
    *pweight_(block) += weight;
    int i = ic3[0] + 2*(ic3[1]+2*ic3[2]);
    *pweight_child_(block,i) = weight;
  }
  if ((!block->is_leaf()) && psync_weight_(block)->next()) {
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
  {
    char buffer[80];
    sprintf (buffer,"recv_index %d %d\n",index,count);
    TRACE_ORDER_BLOCK(buffer,block);
  }
  if (!self) {
    const int rank = cello::rank();
    int na3[3];
    cello::simulation()->hierarchy()->root_blocks(na3,na3+1,na3+2);
    Index index_next = block->index().next(rank,na3,block->is_leaf(),min_level_);
    *pindex_(block) = index;
    *pcount_(block) = count;
    *pnext_(block) = index_next;
  }
  if (psync_index_(block)->next()) {
    {
      char buffer[80];
      sprintf (buffer,"complete %d %d\n",index,count);
      TRACE_ORDER_BLOCK(buffer,block);
    } 
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

long long * MethodOrderMorton::pindex_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_index_);
}

//----------------------------------------------------------------------

long long * MethodOrderMorton::pcount_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_count_);
}

//----------------------------------------------------------------------

Index * MethodOrderMorton::pnext_(Block * block)
{
  Scalar<Index> scalar(cello::scalar_descr_index(),
                     block->data()->scalar_data_index());
  return scalar.value(is_next_);
}

//----------------------------------------------------------------------

Sync * MethodOrderMorton::psync_index_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return scalar.value(is_sync_index_);
}

//----------------------------------------------------------------------

long long * MethodOrderMorton::pweight_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_weight_);
}

//----------------------------------------------------------------------

long long * MethodOrderMorton::pweight_child_(Block * block, int i)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_weight_child_)+i;
}

//----------------------------------------------------------------------

Sync * MethodOrderMorton::psync_weight_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return scalar.value(is_sync_weight_);
}

