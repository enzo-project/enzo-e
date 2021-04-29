// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgOutput.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgOutput Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_MSG_OUTPUT

//----------------------------------------------------------------------

long MsgOutput::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

MsgOutput::MsgOutput()
  : CMessage_MsgOutput(),
    is_local_(true),
    index_send_(),
    block_trace_(),
    method_output_(nullptr),
    data_msg_(nullptr),
    buffer_(nullptr)
{
  ++counter[cello::index_static()]; 
}

//----------------------------------------------------------------------

MsgOutput::MsgOutput (BlockTrace block_trace,MethodOutput * method_output) 
  : CMessage_MsgOutput(),
    is_local_(true),
    index_send_(),
    block_trace_(block_trace),
    method_output_(method_output),
    data_msg_(nullptr),
    buffer_(nullptr)
{
  ++counter[cello::index_static()]; 
}

//----------------------------------------------------------------------

MsgOutput::~MsgOutput()
{
  --counter[cello::index_static()];
  delete data_msg_;
  data_msg_ = 0;
  CkFreeMsg (buffer_);
  buffer_=nullptr;
}

//----------------------------------------------------------------------

void MsgOutput::set_data_msg  (DataMsg * data_msg) 
{
  if (data_msg_) {
    WARNING ("MsgOutput::set_data_msg()",
	     "overwriting existing data_msg_");
    delete data_msg_;
  }
  data_msg_ = data_msg;
}

//----------------------------------------------------------------------

void * MsgOutput::pack (MsgOutput * msg)
{
  // Return with buffer if already packed
  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = 0;

  // determine buffer size

  size += msg->block_trace_.data_size();
  size += sizeof(void *); // method_output_
  
  int have_data = (msg->data_msg_ != nullptr);
  if (have_data) {
    size += msg->data_msg_->data_size();
  }

  // allocate buffer using CkAllocBuffer()

  char * buffer = (char *) CkAllocBuffer (msg,size);

  // serialize message data into buffer 

  union {
    char * pc;
    int  * pi;
    MethodOutput ** pp;
  };

  pc = buffer;

  pc = msg->block_trace_.save_data(pc);
  (*pp++) = msg->method_output_;

  // data_msg_;
  have_data = (msg->data_msg_ != nullptr);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }

  ASSERT2("MsgOutput::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  CkFreeMsg (msg);

  // Return the buffer
  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgOutput * MsgOutput::unpack(void * buffer)
{

  // Allocate storage using CkAllocBuffer (not new!)

  MsgOutput * msg = (MsgOutput *) CkAllocBuffer (buffer,sizeof(MsgOutput));

  msg = new ((void*)msg) MsgOutput;
  
  msg->is_local_ = false;

  // de-serialize message data from input buffer into allocated message

  union {
    char   * pc;
    int    * pi;
    MethodOutput ** pp;
  };

  pc = (char *) buffer;

  // BlockTrace

  pc = msg->block_trace_.load_data(pc);
  msg->method_output_ = (*pp++);
  
  // data_msg_
  int have_data = (*pi++);
  if (have_data) {
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = nullptr;
  }

  // Save the input buffer for freeing later

  msg->buffer_ = buffer;
  
  return msg;
}

//----------------------------------------------------------------------

void MsgOutput::update (Data * data)
{
  // return if no data to update
  if (data_msg_ == nullptr) return;

  data_msg_->update(data,is_local_);

  if (!is_local_) {
      CkFreeMsg (buffer_);
      buffer_ = nullptr;
  }
}

Index MsgOutput::index_send()
{ return index_send_; }

void MsgOutput::set_index_send(Index index)
{ index_send_ = index; }
