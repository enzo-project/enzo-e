// See LICENSE_CELLO file for license and copyright information

/// @file     charm_EnzoMsgCheck.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-19
/// @brief    [\ref Charm] Declaration of the EnzoMsgCheck Charm++ message

#include "enzo.hpp"
#include "charm_enzo.hpp"

#include "Enzo/io/io.hpp" // IoEnzoBlock

//----------------------------------------------------------------------

long EnzoMsgCheck::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

EnzoMsgCheck::EnzoMsgCheck()
  : CMessage_EnzoMsgCheck(),
    is_local_(true),
    index_send_(),
    data_msg_(nullptr),
    buffer_(nullptr),
    block_name_(),
    block_level_(),
    block_lower_(),
    block_upper_(),
    block_size_(),
    tag_(),
    io_block_(),
    index_this_(),
    index_next_(),
    name_this_(),
    name_next_(),
    index_block_(),
    is_first_(),
    is_last_(),
    name_dir_(),
    index_file_(-1),
    order_index_(-1),
    order_count_(-1)
{
  ++counter[cello::index_static()];
  cello::hex_string(tag_,TAG_LEN);
  std::fill_n(adapt_buffer_,ADAPT_BUFFER_SIZE,0.0);
}

//----------------------------------------------------------------------

EnzoMsgCheck::~EnzoMsgCheck()
{
  --counter[cello::index_static()];
  delete data_msg_;
  data_msg_ = nullptr;
  delete io_block_;
  io_block_ = nullptr;
  CkFreeMsg (buffer_);
  buffer_=nullptr;
}

//----------------------------------------------------------------------

void EnzoMsgCheck::set_data_msg  (DataMsg * data_msg)
{
  if (data_msg_) {
    WARNING ("EnzoMsgCheck::set_data_msg()",
	     "overwriting existing data_msg_");
    delete data_msg_;
  }
  data_msg_ = data_msg;
}

//----------------------------------------------------------------------

void * EnzoMsgCheck::pack (EnzoMsgCheck * msg)
{
  // Return with buffer if already packed
  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = msg->size_();

  // determine buffer size

  //--------------------------------------------------

  // allocate buffer using CkAllocBuffer()

  char * buffer = (char *) CkAllocBuffer (msg,size);

  // serialize message data into buffer

  char * pc = buffer;

  pc = msg->save_ (pc);

  ASSERT2("EnzoMsgCheck::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  delete msg;

  // Return the buffer
  return (void *) buffer;
}

//----------------------------------------------------------------------

EnzoMsgCheck * EnzoMsgCheck::unpack(void * buffer)
{

  // Allocate storage using CkAllocBuffer (not new!)

  EnzoMsgCheck * msg = (EnzoMsgCheck *) CkAllocBuffer (buffer,sizeof(EnzoMsgCheck));

  msg = new ((void*)msg) EnzoMsgCheck;

  msg->is_local_ = false;

  // de-serialize message data from input buffer into allocated message

  char * pc = (char *) buffer;

  pc = msg->load_(pc);

  // Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

int EnzoMsgCheck::size_()
{
  int size = 0;
  SIZE_OBJECT_TYPE(size,index_send_);
  SIZE_OBJECT_PTR_TYPE(size,DataMsg,data_msg_);
  SIZE_STRING_TYPE(size,block_name_);
  SIZE_SCALAR_TYPE(size,int,block_level_);
  SIZE_ARRAY_TYPE(size,double,block_lower_,3);
  SIZE_ARRAY_TYPE(size,double,block_upper_,3);
  SIZE_ARRAY_TYPE(size,int,block_size_,3);
  SIZE_ARRAY_TYPE(size,char,tag_,TAG_LEN+1);
  SIZE_OBJECT_PTR_TYPE(size,IoEnzoBlock,io_block_);
  SIZE_OBJECT_TYPE(size,index_this_);
  SIZE_OBJECT_TYPE(size,index_next_);
  SIZE_STRING_TYPE(size,name_this_);
  SIZE_STRING_TYPE(size,name_next_);
  SIZE_SCALAR_TYPE(size,long long,index_block_);
  SIZE_SCALAR_TYPE(size,bool,is_first_);
  SIZE_SCALAR_TYPE(size,bool,is_last_);
  SIZE_STRING_TYPE(size,name_dir_);
  SIZE_SCALAR_TYPE(size,int,index_file_);
  SIZE_ARRAY_TYPE (size,int,adapt_buffer_,ADAPT_BUFFER_SIZE);
  SIZE_SCALAR_TYPE(size,int,order_index_);
  SIZE_SCALAR_TYPE(size,int,order_count_);
  return size;
}

//----------------------------------------------------------------------

char * EnzoMsgCheck::save_(char * pc)
{
  SAVE_OBJECT_TYPE(pc,index_send_);
  SAVE_OBJECT_PTR_TYPE(pc,DataMsg,data_msg_);
  SAVE_STRING_TYPE(pc,block_name_);
  SAVE_SCALAR_TYPE(pc,int,block_level_);
  SAVE_ARRAY_TYPE(pc,double,block_lower_,3);
  SAVE_ARRAY_TYPE(pc,double,block_upper_,3);
  SAVE_ARRAY_TYPE(pc,int,block_size_,3);
  SAVE_ARRAY_TYPE(pc,char,tag_,TAG_LEN+1);
  SAVE_OBJECT_PTR_TYPE(pc,IoEnzoBlock,io_block_);
  SAVE_OBJECT_TYPE(pc,index_this_);
  SAVE_OBJECT_TYPE(pc,index_next_);
  SAVE_STRING_TYPE(pc,name_this_);
  SAVE_STRING_TYPE(pc,name_next_);
  SAVE_SCALAR_TYPE(pc,long long,index_block_);
  SAVE_SCALAR_TYPE(pc,bool,is_first_);
  SAVE_SCALAR_TYPE(pc,bool,is_last_);
  SAVE_STRING_TYPE(pc,name_dir_);
  SAVE_SCALAR_TYPE(pc,int,index_file_);
  SAVE_ARRAY_TYPE (pc,int,adapt_buffer_,ADAPT_BUFFER_SIZE);
  SAVE_SCALAR_TYPE(pc,int,order_index_);
  SAVE_SCALAR_TYPE(pc,int,order_count_);
  return pc;
}

//----------------------------------------------------------------------

char * EnzoMsgCheck::load_(char * pc)
{
  LOAD_OBJECT_TYPE(pc,index_send_);
  LOAD_OBJECT_PTR_TYPE(pc,DataMsg,data_msg_);
  LOAD_STRING_TYPE(pc,block_name_);
  LOAD_SCALAR_TYPE(pc,int,block_level_);
  LOAD_ARRAY_TYPE(pc,double,block_lower_,3);
  LOAD_ARRAY_TYPE(pc,double,block_upper_,3);
  LOAD_ARRAY_TYPE(pc,int,block_size_,3);
  LOAD_ARRAY_TYPE(pc,char,tag_,TAG_LEN+1);
  LOAD_OBJECT_PTR_TYPE(pc,IoEnzoBlock,io_block_);
  LOAD_OBJECT_TYPE(pc,index_this_);
  LOAD_OBJECT_TYPE(pc,index_next_);
  LOAD_STRING_TYPE(pc,name_this_);
  LOAD_STRING_TYPE(pc,name_next_);
  LOAD_SCALAR_TYPE(pc,long long,index_block_);
  LOAD_SCALAR_TYPE(pc,bool,is_first_);
  LOAD_SCALAR_TYPE(pc,bool,is_last_);
  LOAD_STRING_TYPE(pc,name_dir_);
  LOAD_SCALAR_TYPE(pc,int,index_file_);
  LOAD_ARRAY_TYPE (pc,int,adapt_buffer_,ADAPT_BUFFER_SIZE);
  LOAD_SCALAR_TYPE(pc,int,order_index_);
  LOAD_SCALAR_TYPE(pc,int,order_count_);
  return pc;
}
//----------------------------------------------------------------------

void EnzoMsgCheck::update (EnzoBlock * block)
{
  if (io_block_ != nullptr) {
    io_block_->save_to(block);
  }
  update(block->data());
}

//----------------------------------------------------------------------

void EnzoMsgCheck::update (Data * data)
{
  // return if no data to update
  if (data_msg_ == nullptr) return;
  const bool is_kept = false;
  data_msg_->update(data,is_local_,is_kept);

  if (!is_local_) {
    CkFreeMsg (buffer_);
    buffer_ = nullptr;
  }
}

//----------------------------------------------------------------------

Index EnzoMsgCheck::index_send()
{ return index_send_; }

//----------------------------------------------------------------------

void EnzoMsgCheck::set_index_send(Index index)
{ index_send_ = index; }

//----------------------------------------------------------------------

void EnzoMsgCheck::set_block (Block * block)
{
  block_name_ = block->name();
  block_level_ = block->level();
  block->data()->lower(block_lower_,block_lower_+1,block_lower_+2);
  block->data()->upper(block_upper_,block_upper_+1,block_upper_+2);
  block->data()->field().size(block_size_,block_size_+1,block_size_+2);

  io_block_ = (IoEnzoBlock *)enzo::factory()->create_io_block();
  io_block_->set_block(block);
  Scalar<long long> scalar(cello::scalar_descr_long_long(),
                           block->data()->scalar_data_long_long());
  block->get_order(&order_index_, &order_count_);
}

//----------------------------------------------------------------------

void EnzoMsgCheck::del_block()
{
  delete data_msg_;
  data_msg_ = nullptr;
  delete io_block_;
  io_block_ = nullptr;
}

//----------------------------------------------------------------------

void EnzoMsgCheck::print (const char * msg)
{
  CkPrintf ("%d ENZO_MSG_CHECK====================\n",CkMyPe());
  CkPrintf ("%d ENZO_MSG_CHECK tag %s %s\n",CkMyPe(),msg,tag_);
  CkPrintf ("%d ENZO_MSG_CHECK is_local %d\n",CkMyPe(),is_local_);
  int v3[3];
  index_this_.values(v3);
  CkPrintf ("%d ENZO_MSG_CHECK data_msg_ %p\n",CkMyPe(),(void *)data_msg_);
  CkPrintf ("%d ENZO_MSG_CHECK io_block_ %p\n",CkMyPe(),(void *)io_block_);
  CkPrintf ("%d ENZO_MSG_CHECK block_name_ %s\n",CkMyPe(),block_name_.c_str());
  CkPrintf ("%d ENZO_MSG_CHECK block_size_ %d %d %d\n",CkMyPe(),
            block_size_[0],block_size_[1],block_size_[2]);
  CkPrintf ("%d ENZO_MSG_CHECK block_level_ %d\n",CkMyPe(),block_level_);
  CkPrintf ("%d ENZO_MSG_CHECK block_lower %g %g %g\n",CkMyPe(),
            block_lower_[0],block_lower_[1],block_lower_[2]);
  CkPrintf ("%d ENZO_MSG_CHECK block_upper %g %g %g\n",CkMyPe(),
            block_upper_[0],block_upper_[1],block_upper_[2]);
  CkPrintf ("%d ENZO_MSG_CHECK\n",CkMyPe());
  if (io_block_ != nullptr) {
    io_block_->print(msg);
  }
  fflush(stdout);
}
