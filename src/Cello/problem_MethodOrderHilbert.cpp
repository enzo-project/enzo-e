// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOrderHilbert.cpp
/// @author   John Brennan (john.brennnan@mu.ie)
/// @date     2023-11-13
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

MethodOrderHilbert::MethodOrderHilbert(int min_level) throw ()
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

void MethodOrderHilbert::compute (Block * block) throw()
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

  CkCallback callback (CkIndex_Block::r_method_order_hilbert_continue(nullptr),
                       block->proxy_array());

  block->contribute (callback);

}

//----------------------------------------------------------------------

void Block::r_method_order_hilbert_continue(CkReductionMsg * msg)
{
  delete msg;
  static_cast<MethodOrderHilbert*>
    (this->method())->compute_continue(this);
}

//----------------------------------------------------------------------

void MethodOrderHilbert::compute_continue(Block * block)
{
  TRACE_ORDER_BLOCK("continue",block);
  send_weight(block, 0, true);
}

//======================================================================
void MethodOrderHilbert::send_weight(Block * block, int weight_child, bool self)
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
    cello::block_array()[index_parent].p_method_order_hilbert_weight
      (ic3,weight,block->index());
    send_index(block, 0, 0, self);
  } else if (level == min_level_) {

    const int rank = cello::rank();
    int na3[3];
    cello::simulation()->hierarchy()->root_blocks(na3,na3+1,na3+2);

    // TODO: use a "Hilbert_next" function here.
    // Index index_next = block->index().next(rank,na3,block->is_leaf(),min_level_);
    Index index_next = hilbert_next(block->index(), rank, block->is_leaf(), min_level_);

    *pindex_(block) = 0;
    *pcount_(block) = 0;
    *pnext_(block) = index_next;

    send_index(block, 0, weight, self);
    if (!self) {
      CkCallback callback
        (CkIndex_Block::r_method_order_hilbert_complete (nullptr),
         block->proxy_array());
      block->contribute (callback);
    }
  }
}

//----------------------------------------------------------------------

void Block::p_method_order_hilbert_weight(int ic3[3], int weight, Index index_child)
{
  static_cast<MethodOrderHilbert*>
    (this->method())->recv_weight(this, ic3,weight,false);
}

//----------------------------------------------------------------------

void MethodOrderHilbert::recv_weight
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

void MethodOrderHilbert::send_index
(Block * block, int index_parent, int count, bool self)
{
  *pcount_(block) = count;
  if (!block->is_leaf()) {
    int index = *pindex_(block) + 1;

    // TODO: Insert here code to loop over child blocks in hilbert order (I think this is done??)
    int children[8];
    hilbert_children(block, children);

    for (int i=0; i<cello::num_children(); i++) {
      int ic3[3];
      ic3[0] = (children[i] >> 0) & 1;
      ic3[1] = (children[i] >> 1) & 1;
      ic3[2] = (children[i] >> 2) & 1;
      Index index_child = block->index().index_child(ic3,min_level_);
      cello::block_array()[index_child].p_method_order_hilbert_index(index,count);
      index += *pweight_child_(block, children[i]);
    }
  }
}

void Block::p_method_order_hilbert_index(int index, int count)
{
  static_cast<MethodOrderHilbert*>
    (this->method())->recv_index(this, index, count, false);
}

void MethodOrderHilbert::recv_index
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

    // TODO: use a "Hilbert_next" function here.
    // Index index_next = block->index().next(rank,na3,block->is_leaf(),min_level_);
    Index index_next = hilbert_next(block->index(), rank, block->is_leaf(), min_level_);
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
    CkCallback callback (CkIndex_Block::r_method_order_hilbert_complete(nullptr),
                       block->proxy_array());
    block->contribute (callback);
  }
}

//----------------------------------------------------------------------

void Block::r_method_order_hilbert_complete(CkReductionMsg * msg)
{
  delete msg;
  static_cast<MethodOrderHilbert*>
    (this->method())->compute_complete(this);
}

//----------------------------------------------------------------------

void MethodOrderHilbert::compute_complete(Block * block)
{
  // Update Block's index and count
  block->set_order(*pindex_(block),*pcount_(block));
  block->compute_done();
}


//======================================================================

long long * MethodOrderHilbert::pindex_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_index_);
}

//----------------------------------------------------------------------

long long * MethodOrderHilbert::pcount_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_count_);
}

//----------------------------------------------------------------------

Index * MethodOrderHilbert::pnext_(Block * block)
{
  Scalar<Index> scalar(cello::scalar_descr_index(),
                     block->data()->scalar_data_index());
  return scalar.value(is_next_);
}

//----------------------------------------------------------------------

Sync * MethodOrderHilbert::psync_index_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return scalar.value(is_sync_index_);
}

//----------------------------------------------------------------------

long long * MethodOrderHilbert::pweight_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_weight_);
}

//----------------------------------------------------------------------

long long * MethodOrderHilbert::pweight_child_(Block * block, int i)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return scalar.value(is_weight_child_)+i;
}

//----------------------------------------------------------------------

Sync * MethodOrderHilbert::psync_weight_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return scalar.value(is_sync_weight_);
}


//=========================================================

void MethodOrderHilbert::hilbert_children(Block * block, int* children){
    Index index = block->index();
    int T = next_hilbert_state(index);

    for (int ic=0; ic<cello::num_children(); ic++) {
        int hilbert_ind = PHM[T][ic];
        children[hilbert_ind] = ic;
    }
}

int MethodOrderHilbert::next_hilbert_state(Index index){
    int level = index.level();

    // Get binary inds for block. Here the ind for a particular
    // direction consists of the array bits (which is 10 bits) and 
    // the first l tree bits where l is the level of the block.
    // ie, if block is on level 5 then x = AAAAAAAAAATTTTT.

    // I think inds are already in the form I need. The bits are
    // AAAAAAAAAATTTTTTTTTTTTTTTTTTTTLL where the array bits are
    // right-justified and the tree bits are left justified. I just
    // need to ignore the last two sign bits.

    // Array bits for blocks with a negative level are already 
    // appropriately padded with zeros on the left when they're
    // created.

    // All of this should mean we can take the index bits 'as-is'
    // and work with them to get the hilbert indices etc.
    int x = index[0];
    int z = index[0];
    int y = index[0];

    int m = 10 + level;
    int T = 0;
    int TOTAL_BITS = 32;
    int xi, yi, zi, xyzi;

    for (int i = 0; i <= m; i++){
        xi = (x >> (TOTAL_BITS-i)) & 1;
        yi = (y >> (TOTAL_BITS-i)) & 1;
        zi = (z >> (TOTAL_BITS-i)) & 1;

        xyzi = (xi << 2) | (yi << 1) | zi;
        T = PNM[T][xyzi];
    }

    return T;
}

int MethodOrderHilbert::HPM[12][8] = {{0, 1, 3, 2, 6, 7, 5, 4},
                                      {6, 7, 5, 4, 0, 1, 3, 2},
                                      {5, 7, 6, 4, 0, 2, 3, 1},
                                      {3, 1, 0, 2, 6, 4, 5, 7},
                                      {0, 1, 5, 4, 6, 7, 3, 2},
                                      {6, 2, 3, 7, 5, 1, 0, 4},
                                      {5, 1, 0, 4, 6, 2, 3, 7},
                                      {3, 2, 6, 7, 5, 4, 0, 1},
                                      {0, 4, 6, 2, 3, 7, 5, 1},
                                      {6, 4, 0, 2, 3, 1, 5, 7},
                                      {5, 7, 3, 1, 0, 2, 6, 4},
                                      {3, 7, 5, 1, 0, 4, 6, 2}};

int MethodOrderHilbert::HNM[12][8] = {{ 8,  4,  4,  3,  3,  5,  5, 10},
                                      { 9,  5,  5,  2,  2,  4,  4, 11},
                                      { 6, 10, 10,  1,  1,  8,  8,  7},
                                      { 7, 11, 11,  0,  0,  9,  9,  6},
                                      { 8,  0,  0,  6,  6,  1,  1, 11},
                                      { 1,  9,  9,  7,  7, 10, 10,  0},
                                      { 2, 10, 10,  4,  4,  9,  9,  3},
                                      {11,  3,  3,  5,  5,  2,  2,  8},
                                      { 0,  4,  4,  9,  9,  7,  7,  2},
                                      { 5,  1,  1,  8,  8,  3,  3,  6},
                                      { 6,  2,  2, 11, 11,  0,  0,  5},
                                      { 3,  7,  7, 10, 10,  4,  4,  1}};

int MethodOrderHilbert::PHM[12][8] = {{0, 1, 3, 2, 7, 6, 4, 5},
                                      {4, 5, 7, 6, 3, 2, 0, 1},
                                      {4, 7, 5, 6, 3, 0, 2, 1},
                                      {2, 1, 3, 0, 5, 6, 4, 7},
                                      {0, 1, 7, 6, 3, 2, 4, 5},
                                      {6, 5, 1, 2, 7, 4, 0, 3},
                                      {2, 1, 5, 6, 3, 0, 4, 7},
                                      {6, 7, 1, 0, 5, 4, 2, 3},
                                      {0, 7, 3, 4, 1, 6, 2, 5},
                                      {2, 5, 3, 4, 1, 6, 0, 7},
                                      {4, 3, 5, 2, 7, 0, 6, 1},
                                      {4, 3, 7, 0, 5, 2, 6, 1}};

int MethodOrderHilbert::PNM[12][8] = {{ 8,  4,  3,  4, 10,  5,  3,  5},
                                      { 2,  4, 11,  4,  2,  5,  9,  5},
                                      { 1,  7,  8,  8,  1,  6, 10, 10},
                                      {11, 11,  0,  7,  9,  9,  0,  6},
                                      { 8,  0, 11,  1,  6,  0,  6,  1},
                                      {10, 10,  9,  9,  0,  7,  1,  7},
                                      {10, 10,  9,  9,  4,  2,  4,  3},
                                      { 2,  8,  3, 11,  2,  5,  3,  5},
                                      { 0,  2,  9,  9,  4,  7,  4,  7},
                                      { 1,  3,  8,  8,  1,  3,  5,  6},
                                      {11, 11,  0,  2,  5,  6,  0,  2},
                                      {10, 10,  1,  3,  4,  7,  4,  7}};



//##########################################################
Index MethodOrderHilbert::hilbert_next (Index index, int rank, bool is_leaf, int min_level)
{
    // make a copy
    Index index_next = index;

    int level = index.level();
    int ARRAY_BITS = 10;
    int m = ARRAY_BITS + level;
    int states[m+1];
    hilbert_states(index, m, states);
    int T = states[m];

    if (! is_leaf) {
        int ic = HPM[T][0];
        index_next.push_child((ic >> 2) & 1, (ic >> 1) & 1, ic & 1, min_level);

    } else {
        // find the first ancestor (including self) that has child[k]=0 for some k.
        int ic3[3] = {0,0,0};
        int level = index_next.level();
        int xyz = 0;
        if (level > min_level) {
            index_next.child(level, ic3, ic3+1, ic3+2, min_level); // NOTE: sets ic3 to child index of this index in parent.
            xyz = (ic3[0] << 2) + (ic3[1] << 1) + (ic3[0]);
        }

        bool last = (level == min_level) || (PHM[T][xyz] == 7);

        // NOTE: this loop walks up the tree until it finds an index where last isn't true.
        while (level > min_level && last) {
            // NOTE: this repeats the same logic as above but for the parent index.
            index_next = index_next.index_parent(min_level);
            level = index_next.level();
            ic3[0] = ic3[1] = ic3[2] =0;
            if (level>min_level) {
                index_next.child(level, ic3, ic3+1, ic3+2, min_level);
                xyz = (ic3[0] << 2) + (ic3[1] << 1) + (ic3[0]);
            }

            T = states[ARRAY_BITS + level];
            bool last = (level == min_level) || (PHM[T][xyz] == 7);
        }

        if (level == min_level) {
            xyz = (ic3[0] << 2) + (ic3[1] << 1) + (ic3[0]);
            int A_k = PHM[T][xyz];
            A_k = (A_k + 1) % 8;
            xyz = HPM[T][A_k];
        } else {
            xyz = 0;
        }
    
        ic3[0] = (xyz >> 2) & 1;
        ic3[1] = (xyz >> 1) & 1;
        ic3[2] = xyz & 1;
        index_next.set_child(level, ic3[0], ic3[1], ic3[2], min_level);
    }
    return index_next;
}

void MethodOrderHilbert::hilbert_states(Index index, int m, int* states) {
    int TOTAL_BITS = 32;

    // int level = index.level();
    // int ARRAY_BITS = 10;
    // int m = ARRAY_BITS + level;
    // int states[m+1];

    int x = index[0];
    int y = index[1];
    int z = index[2];
    int xi, yi, zi, xyzi, T = 0;

    for (int i = 0; i < m; i++) {
        xi = (x >> (TOTAL_BITS-i)) & 1;
        yi = (y >> (TOTAL_BITS-i)) & 1;
        zi = (z >> (TOTAL_BITS-i)) & 1;
        xyzi = (xi << 2) || (yi << 1) || zi;

        states[i] = T;
        T = PNM[T][xyzi];
    }
    states[m] = T;
}