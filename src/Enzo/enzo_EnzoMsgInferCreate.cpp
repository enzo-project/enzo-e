// See LICENSE_CELLO file for license and copyright information

/// @file     charm_EnzoMsgInferCreate.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-19
/// @brief    [\ref Charm] Declaration of the EnzoMsgInferCreate Charm++ message

#include "enzo.hpp"
#include "charm_enzo.hpp"

//----------------------------------------------------------------------

long EnzoMsgInferCreate::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

EnzoMsgInferCreate::EnzoMsgInferCreate()
  : CMessage_EnzoMsgInferCreate(),
    is_local_(true),
    buffer_(nullptr),
    tag_()
{
  ++counter[cello::index_static()];
  cello::hex_string(tag_,TAG_LEN);
}

//----------------------------------------------------------------------

EnzoMsgInferCreate::~EnzoMsgInferCreate()
{
  --counter[cello::index_static()];
  if (is_local_) {
  }
  CkFreeMsg (buffer_);
  buffer_=nullptr;
}

//----------------------------------------------------------------------

void * EnzoMsgInferCreate::pack (EnzoMsgInferCreate * msg)
{
  // Return with buffer if already packed
  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = 0;

  // determine buffer size

  SIZE_ARRAY_TYPE(size,char,msg->tag_,TAG_LEN+1);

  //--------------------------------------------------

  // allocate buffer using CkAllocBuffer()

  char * buffer = (char *) CkAllocBuffer (msg,size);

  // serialize message data into buffer

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  SAVE_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);

  ASSERT2("EnzoMsgInferCreate::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  CkFreeMsg (msg);
  // Return the buffer
  return (void *) buffer;
}

//----------------------------------------------------------------------

EnzoMsgInferCreate * EnzoMsgInferCreate::unpack(void * buffer)
{

  // Allocate storage using CkAllocBuffer (not new!)

  EnzoMsgInferCreate * msg = (EnzoMsgInferCreate *) CkAllocBuffer (buffer,sizeof(EnzoMsgInferCreate));

  msg = new ((void*)msg) EnzoMsgInferCreate;

  msg->is_local_ = false;

  // de-serialize message data from input buffer into allocated message

  union {
    char   * pc;
    int    * pi;
  };

  pc = (char *) buffer;

  LOAD_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);

  // Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void EnzoMsgInferCreate::print (const char * msg)
{
  CkPrintf ("%d ENZO_MSG_INFER_CREATE====================\n",CkMyPe());
  CkPrintf ("%d ENZO_MSG_INFER_CREATE tag %s %s\n",CkMyPe(),msg,tag_);
  CkPrintf ("%d ENZO_MSG_INFER_CREATE is_local %d\n",CkMyPe(),is_local_);
  fflush(stdout);
    
}
