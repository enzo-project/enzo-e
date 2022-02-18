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
    tag_(),
    io_block_(),
    block_name_(),
    block_lower_(),
    block_upper_()
{
  ++counter[cello::index_static()];
  cello::hex_string(tag_,TAG_LEN);
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
    tag_(),
    io_block_(),
    block_name_(),
    block_lower_(),
    block_upper_()
{
  ++counter[cello::index_static()]; 
  cello::hex_string(tag_,TAG_LEN);
}

//----------------------------------------------------------------------

MsgOutput::~MsgOutput()
{
  --counter[cello::index_static()];
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
  SIZE_STRING_TYPE(size,msg->block_name_);
  SIZE_ARRAY_TYPE(size,double,msg->block_lower_,3);
  SIZE_ARRAY_TYPE(size,double,msg->block_upper_,3);
  
  SIZE_ARRAY_TYPE(size,char,msg->tag_,TAG_LEN+1);

  int have_io = (msg->io_block_ != nullptr);
  SIZE_SCALAR_TYPE(size,int,have_io);
  if (have_io) {
    SIZE_OBJECT_TYPE(size,*(msg->io_block_));
  }

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
  SAVE_STRING_TYPE(pc,msg->block_name_);
  SAVE_ARRAY_TYPE(pc,double,msg->block_lower_,3);
  SAVE_ARRAY_TYPE(pc,double,msg->block_upper_,3);

  SAVE_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);

  have_io = (msg->io_block_ != nullptr);
  SAVE_SCALAR_TYPE(pc,int,have_io);
  
  if (have_io) {
    SAVE_OBJECT_TYPE(pc,*(msg->io_block_));
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
  LOAD_STRING_TYPE(pc,msg->block_name_);
  LOAD_ARRAY_TYPE(pc,double,msg->block_lower_,3);
  LOAD_ARRAY_TYPE(pc,double,msg->block_upper_,3);

  LOAD_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);
  
  int have_io;
  LOAD_SCALAR_TYPE(pc,int,have_io);
  if (have_io) {

    // create the correct IoBlock (IoBlock or IoEnzoBlock
    msg->io_block_ = msg->method_output_->factory()->create_io_block();

    LOAD_OBJECT_TYPE(pc,*(msg->io_block_));
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

//----------------------------------------------------------------------

Index MsgOutput::index_send()
{ return index_send_; }

//----------------------------------------------------------------------

void MsgOutput::set_index_send(Index index)
{ index_send_ = index; }

//----------------------------------------------------------------------

void MsgOutput::set_block (Block * block, const Factory * factory)
{
  block_name_ = block->name();
  block->data()->lower(block_lower_,block_lower_+1,block_lower_+2);
  block->data()->upper(block_upper_,block_upper_+1,block_upper_+2);
  delete io_block_;
  io_block_ = factory->create_io_block();
  io_block_->set_block(block);
}

//----------------------------------------------------------------------

void MsgOutput::del_block()
{
  delete data_msg_;
  data_msg_ = nullptr;
  delete io_block_;
  io_block_ = nullptr;
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
  block_trace_.print(msg);
  CkPrintf ("MSG_OUTPUT method_output_ %p\n",(void *)method_output_);
  CkPrintf ("MSG_OUTPUT file_ %p\n",(void *)file_);
  CkPrintf ("MSG_OUTPUT data_msg_ %p\n",(void *)data_msg_);
  CkPrintf ("MSG_OUTPUT\n");
}
