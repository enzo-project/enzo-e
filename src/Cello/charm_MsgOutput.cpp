// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgOutput.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgOutput Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

//----------------------------------------------------------------------

long MsgOutput::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

MsgOutput::MsgOutput()
  : CMessage_MsgOutput(),
    is_local_(true),
    index_send_(),
    block_trace_(),
    method_output_(nullptr),
    file_(nullptr),
    data_msg_(nullptr),
    buffer_(nullptr),
    block_name(),
    tag_()
{
  ++counter[cello::index_static()];
  cello::hex_string(tag_,TAG_LEN-1);
}

//----------------------------------------------------------------------

MsgOutput::MsgOutput
(
 BlockTrace block_trace,
 MethodOutput * method_output,
 FileHdf5 * file) 
  : CMessage_MsgOutput(),
    is_local_(true),
    index_send_(),
    block_trace_(block_trace),
    method_output_(method_output),
    file_(file),
    data_msg_(nullptr),
    buffer_(nullptr),
    block_name()
{
  cello::hex_string(tag_,TAG_LEN-1);
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

  SIZE_OBJECT_TYPE(size,msg->index_send_);
  SIZE_OBJECT_TYPE(size,msg->block_trace_);
  size += sizeof(void *); // method_output_
  size += sizeof(void *); // file_
  
  // data_msg_;
  int have_data = (msg->data_msg_ != nullptr);
  size += sizeof(int);
  if (have_data) {
    size += msg->data_msg_->data_size();
  }

  // Block name
  SIZE_STRING_TYPE(size,msg->block_name);
  SIZE_ARRAY_TYPE(size,char,msg->tag_,TAG_LEN);

  //--------------------------------------------------
  
  // allocate buffer using CkAllocBuffer()

  char * buffer = (char *) CkAllocBuffer (msg,size);

  // serialize message data into buffer 

  union {
    char * pc;
    int  * pi;
    MethodOutput ** pm;
    FileHdf5 ** pf;
  };

  pc = buffer;

  SAVE_OBJECT_TYPE(pc,msg->index_send_);
  SAVE_OBJECT_TYPE(pc,msg->block_trace_);
  (*pm++) = msg->method_output_;
  (*pf++) = msg->file_;

  // data_msg_;
  have_data = (msg->data_msg_ != nullptr);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }

  // Block name
  SAVE_STRING_TYPE(pc,msg->block_name);
  SAVE_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN);

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
    MethodOutput ** pm;
    FileHdf5 ** pf;
  };

  pc = (char *) buffer;

  LOAD_OBJECT_TYPE(pc,msg->index_send_);
  LOAD_OBJECT_TYPE(pc,msg->block_trace_);
  msg->method_output_ = (*pm++);
  msg->file_          = (*pf++);
  
  // data_msg_
  int have_data = (*pi++);
  if (have_data) {
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = nullptr;
  }

  // Block name
  LOAD_STRING_TYPE(pc,msg->block_name);
  LOAD_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN);
  
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

//----------------------------------------------------------------------

Index MsgOutput::index_send()
{ return index_send_; }

//----------------------------------------------------------------------

void MsgOutput::set_index_send(Index index)
{ index_send_ = index; }

//----------------------------------------------------------------------

void MsgOutput::set_block (Block * block)
{
  block_name = block->name();
}

//----------------------------------------------------------------------

void MsgOutput::del_block()
{
}

//----------------------------------------------------------------------

void MsgOutput::print (const char * msg)
{
    CkPrintf ("MSG_OUTPUT====================\n");
    CkPrintf ("MSG_OUTPUT tag %s %s\n",msg,tag_);
    CkPrintf ("MSG_OUTPUT is_local %d\n",is_local_);
    int v3[3];
    index_send_.values(v3);
    CkPrintf ("MSG_OUTPUT index_send values %d %d %d\n",v3[0],v3[1],v3[2]);
    CkPrintf ("MSG_OUTPUT block_trace_ %p\n",block_trace_);
    block_trace_.print(msg);
    CkPrintf ("MSG_OUTPUT method_output_ %p\n",method_output_);
    CkPrintf ("MSG_OUTPUT file_ %p\n",file_);
    CkPrintf ("MSG_OUTPUT data_msg_ %p\n",data_msg_);
    CkPrintf ("MSG_OUTPUT\n");
  }
