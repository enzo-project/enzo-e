// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgRefresh Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_NEW_REFRESH
// #define TRACE_NEW_REFRESH

//----------------------------------------------------------------------

int MsgRefresh::id_count = -1;
long MsgRefresh::counter = 0;

//----------------------------------------------------------------------

MsgRefresh::MsgRefresh()
    : CMessage_MsgRefresh(),
      is_local_(true),
      id_(-1),
      data_msg_(NULL),
      buffer_(NULL)
{  ++counter; }

//----------------------------------------------------------------------

MsgRefresh::~MsgRefresh()
{
  --counter;
  // delete data_msg_;
  // data_msg_ = 0;
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
#ifdef TRACE_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgRefresh::pack()\n",msg);
#endif
  // Update ID (for debugging)

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();

  int size = 0;

  size += sizeof(int); // id_

  size += sizeof(int); // have_data

  int have_data = (msg->data_msg_ != NULL);
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

  (*pi++) = msg->id_;

  have_data = (msg->data_msg_ != NULL);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }

  delete msg;

  // Return the buffer

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p MsgRefresh pack message size %d %d\n",msg,(pc - (char*)buffer),size);
#endif
  ASSERT2("MsgRefresh::pack()",
	  "buffer size mismatch %d allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgRefresh * MsgRefresh::unpack(void * buffer)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  MsgRefresh * msg = 
    (MsgRefresh *) CkAllocBuffer (buffer,sizeof(MsgRefresh));

  msg = new ((void*)msg) MsgRefresh;
#ifdef TRACE_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgRefresh::unpack()\n",msg);
#endif
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgRefresh::unpack()\n",msg);
#endif
  
  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  union {
    char * pc;
    int  * pi;
  };

  pc = (char *) buffer;

  msg->id_ = (*pi++);

  int have_data = (*pi++);
  if (have_data) {
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = NULL;
  }

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p MsgRefresh unpack message size %d\n",msg,(pc - (char*)buffer));
#endif

  return msg;
}

//----------------------------------------------------------------------

void MsgRefresh::update (Data * data)
{
  if (data_msg_ == NULL) return;

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgRefresh::update()\n",this);
#endif
#ifdef TRACE_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgRefresh::update()\n",this);
#endif
  data_msg_->update(data,is_local_);
  if (!is_local_) {
      CkFreeMsg (buffer_);
  }
}
