// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgAdapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgAdapt Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

//----------------------------------------------------------------------

long MsgAdapt::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

MsgAdapt::~MsgAdapt()
{
  --counter[cello::index_static()];

  CkFreeMsg (buffer_);
  buffer_=nullptr;
}

//----------------------------------------------------------------------

void * MsgAdapt::pack (MsgAdapt * msg)
{
  if (msg->buffer_ != NULL) return msg->buffer_;

  int size = 0;

  SIZE_SCALAR_TYPE(size,int,   msg->adapt_step_);
  SIZE_SCALAR_TYPE(size,Index, msg->index_);
  SIZE_ARRAY_TYPE (size,int,   msg->ic3_,3);
  SIZE_VECTOR_TYPE(size,int,   msg->ofv_[0]);
  SIZE_VECTOR_TYPE(size,int,   msg->ofv_[1]);
  SIZE_VECTOR_TYPE(size,int,   msg->ofv_[2]);
  SIZE_SCALAR_TYPE(size,int,   msg->level_now_);
  SIZE_SCALAR_TYPE(size,int,   msg->level_min_);
  SIZE_SCALAR_TYPE(size,int,   msg->level_max_);
  SIZE_SCALAR_TYPE(size,bool,  msg->can_coarsen_);
  SIZE_ARRAY_TYPE(size,char,msg->tag_,TAG_LEN+1);

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

  SAVE_SCALAR_TYPE(pc,int,   msg->adapt_step_);
  SAVE_SCALAR_TYPE(pc,Index, msg->index_);
  SAVE_ARRAY_TYPE (pc,int,   msg->ic3_,3);
  SAVE_VECTOR_TYPE(pc,int,   msg->ofv_[0]);
  SAVE_VECTOR_TYPE(pc,int,   msg->ofv_[1]);
  SAVE_VECTOR_TYPE(pc,int,   msg->ofv_[2]);
  SAVE_SCALAR_TYPE(pc,int,   msg->level_now_);
  SAVE_SCALAR_TYPE(pc,int,   msg->level_min_);
  SAVE_SCALAR_TYPE(pc,int,   msg->level_max_);
  SAVE_SCALAR_TYPE(pc,bool,  msg->can_coarsen_);
  SAVE_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);
  
  ASSERT2("MsgRefresh::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  delete msg;

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgAdapt * MsgAdapt::unpack(void * buffer)
{

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  MsgAdapt * msg = 
    (MsgAdapt *) CkAllocBuffer (buffer,sizeof(MsgAdapt));

  msg = new ((void*)msg) MsgAdapt;
  
  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  union {
    char * pc;
    int  * pi;
  };

  pc = (char *) buffer;

  LOAD_SCALAR_TYPE(pc,int,   msg->adapt_step_);
  LOAD_SCALAR_TYPE(pc,Index, msg->index_);
  LOAD_ARRAY_TYPE (pc,int,   msg->ic3_,3);
  LOAD_VECTOR_TYPE(pc,int,   msg->ofv_[0]);
  LOAD_VECTOR_TYPE(pc,int,   msg->ofv_[1]);
  LOAD_VECTOR_TYPE(pc,int,   msg->ofv_[2]);
  LOAD_SCALAR_TYPE(pc,int,   msg->level_now_);
  LOAD_SCALAR_TYPE(pc,int,   msg->level_min_);
  LOAD_SCALAR_TYPE(pc,int,   msg->level_max_);
  LOAD_SCALAR_TYPE(pc,bool,  msg->can_coarsen_);
  LOAD_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void MsgAdapt::print (const char * msg)
{
  const int ip = CkMyPe();
  //  CkPrintf ("%d MSG_ADAPT====================\n",ip);
  CkPrintf ("%d MSG_ADAPT tag %s %p %s\n",ip,msg,(void *)this,tag_);
  //  CkPrintf ("%d MSG_ADAPT adapt_step %d\n",ip,adapt_step_);
  //  CkPrintf ("%d MSG_ADAPT ic3        %d %d %d\n",ip,ic3_[0],ic3_[1],ic3_[2]);
  //  CkPrintf ("%d MSG_ADAPT of3        %d %d %d\n",ip,of3_[0],of3_[1],of3_[2]);
  //  CkPrintf ("%d MSG_ADAPT level_now  %d\n",ip,level_now_);
  //  CkPrintf ("%d MSG_ADAPT level_min  %d\n",ip,level_min_);
  //  CkPrintf ("%d MSG_ADAPT level_max  %d\n",ip,level_max_);
  //  CkPrintf ("%d MSG_ADAPT can_coarsen %d\n",ip,can_coarsen_);
  //  CkPrintf ("%d \n",ip);
  fflush(stdout);
}
