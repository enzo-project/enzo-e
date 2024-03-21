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
  is_sync_index_   = cello::scalar_descr_sync()->new_value(name() + ":sync_index");
  is_sync_weight_  = cello::scalar_descr_sync()->new_value(name() + ":sync_weight");
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

    int children[cello::num_children()];
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
    int ARRAY_BITS = 10;
    int m = ARRAY_BITS + index.level();;
    int states[m+1];
    hilbert_states(index, m, states);
    int T = states[m];

    for (int ic=0; ic<cello::num_children(); ic++) {
        int hilbert_ind = coord_to_hilbert_ind(T, ic);
        children[hilbert_ind] = ic;
    }
}

Index MethodOrderHilbert::hilbert_next (Index index, int rank, bool is_leaf, int min_level)
{
    Index index_next = index;
    int level = index.level();
    int ARRAY_BITS = 10;
    int m = ARRAY_BITS + level;
    int states[m+1];
    hilbert_states(index, m, states);

    // If the block has children then the first child (according to the hilbert order)
    // is the next index.
    if (! is_leaf) {
        int T = states[m];
        int ic = hilbert_ind_to_coord(T, 0);
        index_next.push_child(ic & 1, (ic >> 1) & 1, (ic >> 2) & 1, min_level);

    // Otherwise the next index is the next child of the first ancestor block whose children
    // haven't been completely processed yet
    } else {
        // find the first ancestor (including self) that has child[k]=0 for some k.
        int ic3[3] = {0,0,0};
        int level = index_next.level();
        int zyx = 0;
        if (level > min_level) {
            // NOTE: sets ic3 to child index of this index in parent.
            index_next.child(level, ic3, ic3+1, ic3+2, min_level);
            zyx = (ic3[2] << 2) + (ic3[1] << 1) + (ic3[0]);
        }

        int T = states[ARRAY_BITS + level - 1];
        bool last = (level == min_level) || is_last_child(T, zyx);

        // NOTE: this loop walks up the tree until it finds an index where last isn't true.
        while (level > min_level && last) {

            // NOTE: this repeats the same logic as above but for the parent index.
            index_next = index_next.index_parent(min_level);
            level = index_next.level();
            ic3[0] = ic3[1] = ic3[2] =0;
            if (level>min_level) {
                index_next.child(level, ic3, ic3+1, ic3+2, min_level);
                zyx = (ic3[2] << 2) + (ic3[1] << 1) + (ic3[0]);
            }

            T = states[ARRAY_BITS + level - 1];
            last = (level == min_level) || is_last_child(T, zyx);
        }

        T = states[ARRAY_BITS + level - 1];
        zyx = (ic3[2] << 2) + (ic3[1] << 1) + (ic3[0]);
        int A_k = coord_to_hilbert_ind(T, zyx);
        A_k = (A_k + 1) % 8;
        zyx = hilbert_ind_to_coord(T, A_k);
    
        ic3[2] = (zyx >> 2) & 1;
        ic3[1] = (zyx >> 1) & 1;
        ic3[0] = (zyx >> 0) & 1;
        index_next.set_child(level, ic3[0], ic3[1], ic3[2], min_level);
    }
    return index_next;
}

void MethodOrderHilbert::hilbert_states(Index index, int m, int* states) {
    // Note: the length of the int array 'states' is expected to be m+1.
    int TOTAL_BITS = 32;
    int xi, yi, zi, zyxi, T = 0;

    int array_bits[3] = {0, 0, 0}, tree_bits[3] = {0, 0, 0};
    index.array(array_bits, array_bits+1, array_bits+2);
    index.tree(tree_bits, tree_bits+1, tree_bits+2);
    int shift = index.level() < 0 ? -1 * index.level() : 0;
    int x = ((array_bits[0] << (20 + shift)) | tree_bits[0]) << 2;
    int y = ((array_bits[1] << (20 + shift)) | tree_bits[1]) << 2;
    int z = ((array_bits[2] << (20 + shift)) | tree_bits[2]) << 2;

    for (int i = 0; i < m; i++) {
        xi = (x >> (TOTAL_BITS-i-1)) & 1;
        yi = (y >> (TOTAL_BITS-i-1)) & 1;
        zi = (z >> (TOTAL_BITS-i-1)) & 1;
        zyxi = (zi << 2) | (yi << 1) | xi;

        states[i] = T;
        T = coord_to_next_state(T, zyxi);
    }
    states[m] = T;
}

//========== Hilbert lookup functions ==========

bool MethodOrderHilbert::is_last_child(int T, int coord) {
    int rank = cello::rank();
    int last_ind = 1;
    if (rank == 3) {
      last_ind = 7;
    } else if (rank == 2) {
      last_ind = 3;
    }
    return coord_to_hilbert_ind(T, coord) == last_ind;
}

int MethodOrderHilbert::coord_to_hilbert_ind(int state, int coord) {
    int rank = cello::rank();
    int hilbert_ind;

    if (rank == 3) {
        hilbert_ind = PHM[state][coord & 7];
    } else if (rank == 2) {
        hilbert_ind = CHM[state][coord & 3];
    } else if (rank == 1) {
        hilbert_ind = coord & 1;
    }

    return hilbert_ind;
}

int MethodOrderHilbert::coord_to_next_state(int state, int coord) {
    int rank = cello::rank();
    int next_state;

    if (rank == 3) {
        next_state = PNM[state][coord & 7];
    } else if (rank == 2) {
        next_state = CSM[state][coord & 3];
    } else if (rank == 1) {
        next_state = 0;
    }

    return next_state;
}

int MethodOrderHilbert::hilbert_ind_to_coord(int state, int hilbert_ind) {
    int rank = cello::rank();
    int coord;

    if (rank == 3) {
        coord = HPM[state][hilbert_ind & 7];
    } else if (rank == 2) {
        coord = HCM[state][hilbert_ind & 3];
    } else if (rank == 1) {
        coord = hilbert_ind & 1;
    }

    return coord;
}


//========== 3D Hilbert lookup tables ==========

// Hilbert ind to Point Map
int MethodOrderHilbert::HPM[12][8] = {
    {0, 1, 3, 2, 6, 7, 5, 4},
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

// Hilbert ind to Next state Map
int MethodOrderHilbert::HNM[12][8] = {
    { 8,  4,  4,  3,  3,  5,  5, 10},
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

// Point to Hilbert ind Map
int MethodOrderHilbert::PHM[12][8] = {
    {0, 1, 3, 2, 7, 6, 4, 5},
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

// Point to Next state Map
int MethodOrderHilbert::PNM[12][8] = {
    { 8,  4,  3,  4, 10,  5,  3,  5},
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


//========== 2D Hilbert lookup tables ==========

// Hilbert ind to Coord Map
int MethodOrderHilbert::HCM[4][4] = {
    {0, 1, 3, 2},
    {0, 2, 3, 1},
    {3, 2, 0, 1},
    {3, 1, 0, 2}};

// Hilbert ind to State Map
int MethodOrderHilbert::HSM[4][4] = {
    {1, 0, 0, 3},
    {0, 1, 1, 2},
    {3, 2, 2, 1},
    {2, 3, 3, 0}};

// Coord to Hilbert ind Map
int MethodOrderHilbert::CHM[4][4] = {
    {0, 1, 3, 2},
    {0, 3, 1, 2},
    {2, 3, 1, 0},
    {2, 1, 3, 0}};

// Coord to State Map
int MethodOrderHilbert::CSM[4][4] = {
    {1, 0, 3, 0},
    {0, 2, 1, 1},
    {2, 1, 2, 3},
    {3, 3, 0, 2}};