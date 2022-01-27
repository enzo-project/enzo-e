// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgRefresh Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #undef TRACE_MSG_REFRESH

//----------------------------------------------------------------------

long MsgRefresh::counter[CONFIG_NODE_SIZE] = { };

//----------------------------------------------------------------------

MsgRefresh::MsgRefresh()
    : CMessage_MsgRefresh(),
      is_local_(true),
      id_refresh_(-1),
      data_msg_(nullptr),
      buffer_(nullptr)
#ifdef TRACE_MSG_REFRESH      
    ,
      name_block_(),
      name_type_()
#endif      
{
  ++counter[cello::index_static()];
}

//----------------------------------------------------------------------

MsgRefresh::~MsgRefresh()
{
  --counter[cello::index_static()];
  delete data_msg_;
  data_msg_ = nullptr;
  CkFreeMsg (buffer_);
  buffer_=nullptr;
}

//----------------------------------------------------------------------

void MsgRefresh::set_data_msg  (DataMsg * data_msg) 
{
  if (data_msg_) {
    WARNING ("MsgRefresh::set_data_msg()",
	     "overwriting existing data_msg_");
    delete data_msg_;
  }
  data_msg_ = data_msg;
}

//----------------------------------------------------------------------

void * MsgRefresh::pack (MsgRefresh * msg)
{
  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = 0;

  size += sizeof(int); // id_refresh
  
#ifdef TRACE_MSG_REFRESH      
  size += sizeof(int) + msg->name_block_.size()*sizeof(char);
  size += sizeof(int) + msg->name_type_.size()*sizeof(char);
#endif
  
  size += sizeof(int);  // have_data
  int have_data = (msg->data_msg_ != nullptr);
  
  if (have_data) {
    // data_msg_
    size += msg->data_msg_->data_size();
  }

  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

  char * buffer = (char *) CkAllocBuffer (msg,size);

  //--------------------------------------------------
  //  3. serialize message data into buffer 
  //--------------------------------------------------

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  (*pi++) = msg->id_refresh_;

#ifdef TRACE_MSG_REFRESH
  int n=msg->name_block_.size();
  (*pi++) = n;
  for (int i=0; i<n; i++) (*pc++) = msg->name_block_[i];
  n=msg->name_type_.size();
  (*pi++) = n;
  for (int i=0; i<n; i++) (*pc++) = msg->name_type_[i];
#endif  
  
  have_data = (msg->data_msg_ != nullptr);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }

  delete msg;

  // Return the buffer

  ASSERT2("MsgRefresh::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgRefresh * MsgRefresh::unpack(void * buffer)
{

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  MsgRefresh * msg = 
    (MsgRefresh *) CkAllocBuffer (buffer,sizeof(MsgRefresh));

  msg = new ((void*)msg) MsgRefresh;
  
  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  union {
    char * pc;
    int  * pi;
  };

  pc = (char *) buffer;

  msg->id_refresh_ = (*pi++) ;

#ifdef TRACE_MSG_REFRESH      
  int n = (*pi++);
  msg->name_block_="";
  for (int i=0; i<n; i++) {
    const std::string s(1,(*pc++));
    msg->name_block_.append(s);
  }

  n = (*pi++);
  msg->name_type_="";
  for (int i=0; i<n; i++) {
    const std::string s(1,(*pc++));
    msg->name_type_.append(s);
  }
#endif  
  
  int have_data = (*pi++);
  if (have_data) {
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = nullptr;
  }

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void MsgRefresh::update (Data * data)
{
  if (data_msg_ == nullptr) return;

  data_msg_->update(data,is_local_);

  if (!is_local_) {
      CkFreeMsg (buffer_);
      buffer_ = nullptr;
  }
}

//----------------------------------------------------------------------

void MsgRefresh::print (const char * message)
{
  if (data_msg_) {
    data_msg_->print(message);
  } else {
    CkPrintf ("MSG_REFRESH data_msg_ = nil\n");
  }
}
