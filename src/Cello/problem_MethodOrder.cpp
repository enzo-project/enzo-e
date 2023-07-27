// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOrder.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-01
/// @brief

#include "problem.hpp"

// #define TRACE_ORDER

#ifdef TRACE_ORDER
#  define TRACE_ORDER_BLOCK(MSG,BLOCK)                  \
  CkPrintf ("TRACE_ORDER %d %s %s %lld %lld %g %g\n",      \
            CkMyPe(), std::string(MSG).c_str(),             \
            BLOCK->name().c_str(),                      \
            index_(BLOCK),                              \
            count_(BLOCK),                              \
            windex_(BLOCK),                             \
            wcount_(BLOCK));                            \
  fflush(stdout);
#  define TRACE_BLOCK(MSG,BLOCK)                \
  CkPrintf ("TRACE_ORDER %s %s\n",              \
            std::string(MSG).c_str(),           \
            BLOCK->name().c_str());             \
    fflush(stdout);

#else
#  define TRACE_ORDER_BLOCK(MSG,BLOCK) /* ... */
#  define TRACE_BLOCK(MSG,BLOCK) /* ... */
#endif


//----------------------------------------------------------------------

MethodOrder::MethodOrder(std::string type,
                         int min_level) throw ()
  : Method(),
    is_index_(-1),
    is_count_(-1),
    is_windex_(-1),
    is_wcount_(-1),
    is_next_(-1),
    is_count_child_(-1),
    is_wcount_child_(-1),
    is_sync_index_(-1),
    is_sync_count_(-1),
    min_level_(min_level)
{

  if (type == "morton") {
    type_ = Type::morton;
  } else {
    ERROR1("MethodOrder()",
          "Unsupported MethodOrder type %s",
           type.c_str());
  }
  Refresh * refresh = cello::refresh(ir_post_);
  cello::simulation()->refresh_set_name(ir_post_,name());
  refresh->add_field("density");

  /// Create Scalar data for ordering index
  const int n = cello::num_children();
  auto sd_ll  = cello::scalar_descr_long_long();
  auto sd_ind = cello::scalar_descr_index();
  auto sd_sync = cello::scalar_descr_sync();
  auto sd_d    = cello::scalar_descr_double();
  // output scalars
  is_index_        = sd_ll->  new_value(name() + ":index"); // block index 0..# participating blocks - 1
  is_count_        = sd_ll->  new_value(name() + ":count"); // number of participating blocks
  is_windex_       = sd_d->   new_value(name() + ":windex"); // weighted index
  is_wcount_       = sd_d->   new_value(name() + ":wcount"); // weighted count
  // local scalars
  is_next_         = sd_ind-> new_value(name() + ":next");
  is_count_child_  = sd_ll->  new_value(name() + ":count_child",n);
  is_wcount_child_ = sd_d->   new_value(name() + ":wcount_child",n);
  is_sync_index_   = sd_sync->new_value(name() + ":sync_index");
  is_sync_count_   = sd_sync->new_value(name() + ":sync_count");
}

//======================================================================

void MethodOrder::compute (Block * block) throw()
{
  double weight = weight_(block);
  // call accum_count() and accum_index() on self to initialize sync counters
  int ic3[3];
  block->index().child(block->level(),ic3,ic3+1,ic3+2,min_level_);
  accum_count(block,1,weight,ic3,1 + cello::num_children(block));
  accum_index(block,0,0,0,0,2);
}

//----------------------------------------------------------------------

void MethodOrder::accum_count
(Block * block, int count, double wcount, const int ic3[3], int sync_stop)
{
  // set sync counter when available
  if (sync_stop != 0) {
    sync_count_(block).set_stop(sync_stop);
  } else {
    const int i = child_order_(ic3);
    count_child_(block,i) = count;
    wcount_child_(block,i) = wcount;
  }

  // accumulate subtree block count and weight
  count_(block)  += count;
  wcount_(block) += wcount;

  TRACE_ORDER_BLOCK("accum_count",block);

  // when all counts received, forward total to parent if available, else
  // done counting so switch to accum_index()
  if (sync_count_(block).next()) {
    TRACE_ORDER_BLOCK("accum_count next",block);
    const int level = block->level();
    if (level > min_level_) {
      // create order message and forward to parent
      MsgOrder * msg_order = new MsgOrder;
      int ic3[3];
      block->index().child(level,ic3,ic3+1,ic3+2,min_level_);
      msg_order->set_child(ic3);
      msg_order->set_count(count_(block),wcount_(block));
      auto index_parent = block->index().index_parent(min_level_);
      cello::block_array()[index_parent].p_method_order_accum_count(msg_order);
    } else {
      // else reached coarsest block; switch from accum_count() to
      // accum_index()
      accum_index(block,0,count_(block),0,wcount_(block));
    }
    // reset count sync for next call
    sync_count_(block).reset();
  }
}

//----------------------------------------------------------------------

void Block::p_method_order_accum_count(MsgOrder * msg)
{
  TRACE_BLOCK("p_accum_count",this);
  // unpack order message and delete
  int count;
  double wcount;
  int ic3[3];
  msg->get_count(count,wcount);
  msg->get_child(ic3);
  delete msg;

  // return to accum_count()
  static_cast<MethodOrder*>
    (this->method())->accum_count(this,count, wcount, ic3);
}

//----------------------------------------------------------------------

void MethodOrder::accum_index
(Block * block, int index, int count, double windex, double wcount, int sync_stop)
{
  // set sync counter when available
  if (sync_stop != 0) {
    sync_index_(block).set_stop(sync_stop);
  } else {
    index_(block) = index;
    windex_(block) = windex;
    count_(block) = count;
    wcount_(block) = wcount;
  }

  TRACE_ORDER_BLOCK("accum_index",block);
  // when done, forward to children if available, else accum
  if (sync_index_(block).next()) {
    TRACE_ORDER_BLOCK("accum_index next",block);
    // create order message and forward to children
    int index_order = index_(block);
    double index_worder = windex_(block);
    const int rank = cello::rank();
    const int nc = cello::num_children();
    int ic3[3];
    index += 1;
    windex += weight_(block);
    // loop through children calling accum_index()
    for (int ic=0; ic<nc; ic++) {
      child_order_(ic3,ic);
      auto index_child = block->index().index_child(ic3,min_level_);
      MsgOrder * msg_order = new MsgOrder;
      msg_order->set_index(index, windex);
      msg_order->set_count(count, wcount);
      cello::block_array()[index_child].p_method_order_accum_index(msg_order);
      index  +=  count_child_(block,ic);
      windex += wcount_child_(block,ic);
    }
    // reset index sync for next call
    sync_index_(block).reset();
    compute_complete_(block);
  }
}

//----------------------------------------------------------------------

void Block::p_method_order_accum_index(MsgOrder * msg)
{
  TRACE_BLOCK("p_accum_index",this);
  // unpack order message and delete
  int index, count;
  double windex, wcount;
  msg->get_index(index,windex);
  msg->get_count(count,wcount);
  delete msg;

  // return to accum_index()
  static_cast<MethodOrder*>
    (this->method())->accum_index(this,index, count, windex, wcount);
}

//----------------------------------------------------------------------

void MethodOrder::compute_complete_(Block * block)
{
  TRACE_ORDER_BLOCK("compute_complete",block);
  // Update Block's index and count
  block->set_order(index_(block),count_(block));
  ASSERT2("compute_complete","index %d is not between 0 and count %d\n",index_(block),count_(block),
          0 <= index_(block) && index_(block) < count_(block));
  index_(block) = 0;
  count_(block) = 0;
  windex_(block) = 0.0;
  wcount_(block) = 0.0;
  block->compute_done();
}


//======================================================================

long long & MethodOrder::index_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return *scalar.value(is_index_);
}

//----------------------------------------------------------------------

long long & MethodOrder::count_(Block * block)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return *scalar.value(is_count_);
}

//----------------------------------------------------------------------

double & MethodOrder::windex_(Block * block)
{
  Scalar<double> scalar(cello::scalar_descr_double(),
                     block->data()->scalar_data_double());
  return *scalar.value(is_index_);
}

//----------------------------------------------------------------------

double & MethodOrder::wcount_(Block * block)
{
  Scalar<double> scalar(cello::scalar_descr_double(),
                     block->data()->scalar_data_double());
  return *scalar.value(is_wcount_);
}

//----------------------------------------------------------------------

Index & MethodOrder::next_(Block * block)
{
  Scalar<Index> scalar(cello::scalar_descr_index(),
                     block->data()->scalar_data_index());
  return *scalar.value(is_next_);
}

//----------------------------------------------------------------------

Sync & MethodOrder::sync_index_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return *scalar.value(is_sync_index_);
}

//----------------------------------------------------------------------

long long & MethodOrder::count_child_(Block * block, int index)
{
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                     block->data()->scalar_data_long_long());
  return *(scalar.value(is_count_child_) + index);
}

//----------------------------------------------------------------------

double & MethodOrder::wcount_child_(Block * block, int index)
{
  Scalar<double> scalar(cello::scalar_descr_double(),
                     block->data()->scalar_data_double());
  return *(scalar.value(is_wcount_child_) + index);
}

//----------------------------------------------------------------------

Sync & MethodOrder::sync_count_(Block * block)
{
  Scalar<Sync> scalar(cello::scalar_descr_sync(),
                      block->data()->scalar_data_sync());
  return *scalar.value(is_sync_count_);
}
