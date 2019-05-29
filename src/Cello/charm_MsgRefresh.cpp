// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgRefresh Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_MSG_REFRESH

//----------------------------------------------------------------------

long MsgRefresh::counter[CONFIG_NODE_SIZE] = { };

//----------------------------------------------------------------------

MsgRefresh::MsgRefresh()
    : CMessage_MsgRefresh(),
      is_local_(true),
#ifdef NEW_REFRESH  
      id_refresh_(-1),
#endif      
      data_msg_(nullptr),
      buffer_(nullptr)
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
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("%d %s:%d DEBUG_MSG_REFRESH packing %p\n",
	    CkMyPe(),__FILE__,__LINE__,msg);
#endif  
  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = 0;

#ifdef NEW_REFRESH
  size += sizeof(int); // id_refresh
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

#ifdef NEW_REFRESH
  (*pi++) = msg->id_refresh_;
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("DEBUG_MSG_REFRESH MsgRefresh::pack id_refresh=%d\n",msg->id_refresh_);
#endif  
#endif  

  have_data = (msg->data_msg_ != nullptr);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }

  delete msg;

  // Return the buffer

  ASSERT2("MsgRefresh::pack()",
	  "buffer size mismatch %d allocated %d packed",
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
  
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("%d %s:%d DEBUG_MSG_REFRESH unpacking %p\n",
	    CkMyPe(),__FILE__,__LINE__,msg);
#endif  

  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  union {
    char * pc;
    int  * pi;
  };

  pc = (char *) buffer;

#ifdef NEW_REFRESH
  msg->id_refresh_ = (*pi++) ;
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("DEBUG_MSG_REFRESH MsgRefresh::pack id_refresh=%d\n",msg->id_refresh_);
#endif  
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
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("%d %s:%d DEBUG_MSG_REFRESH updating %p\n",
	    CkMyPe(),__FILE__,__LINE__,this);
#endif  
  if (data_msg_ == nullptr) return;

  data_msg_->update(data,is_local_);

  if (!is_local_) {
      CkFreeMsg (buffer_);
      buffer_ = nullptr;
  }
}
