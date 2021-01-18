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
      id_refresh_(-1),
      data_msg_(nullptr),
      buffer_(nullptr)
#ifdef TRACE_MSG_REFRESH      
    ,
      name_block_(),
      name_type_()
#endif      
{
  cello::hex_string(tag_,8);
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
#ifdef DEBUG_MSG_REFRESH  
  CkPrintf ("DEBUG_CHARM %p set_data_msg %p\n",this,data_msg_);
#endif  
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

  size += 8*sizeof(char);
  
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
  
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("DEBUG_MSG_REFRESH MsgRefresh::pack id_refresh=%d\n",msg->id_refresh_);
#endif  

  have_data = (msg->data_msg_ != nullptr);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }
#ifdef DEBUG_MSG_REFRESH  
  CkPrintf ("DEBUG_CHARM %p pack data_msg_ %p\n",msg,msg->data_msg_);
#endif  

  strncpy(pc,msg->tag_,8);
  msg->tag_[8] = '\0';
  pc += 8*sizeof(char);
  
#ifdef DEBUG_MSG_REFRESH  
  msg->print(std::string(msg->tag_)+":pack");
#endif

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
  
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("DEBUG_MSG_REFRESH MsgRefresh::pack id_refresh=%d\n",msg->id_refresh_);
#endif  

  int have_data = (*pi++);
  if (have_data) {
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = nullptr;
  }
#ifdef DEBUG_MSG_REFRESH  
  CkPrintf ("DEBUG_CHARM %p unpack data_msg_ %p\n",msg,msg->data_msg_);
#endif

  strncpy (msg->tag_,pc,8);
  msg->tag_[8] = '\0';
  pc+=8*sizeof(char);
  
#ifdef DEBUG_MSG_REFRESH  
  CkPrintf ("DEBUG_CHARM %p unpack tag_ %s\n",msg,msg->tag_);
#endif
  
  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

#ifdef DEBUG_MSG_REFRESH  
  msg->print(std::string(msg->tag_)+":unpack");
#endif

  return msg;
}

//----------------------------------------------------------------------

void MsgRefresh::update (Data * data)
{
#ifdef DEBUG_MSG_REFRESH
  CkPrintf ("%d %s:%d DEBUG_MSG_REFRESH updating %p %s\n",
	    CkMyPe(),__FILE__,__LINE__,this,tag_);
#endif  
#ifdef DEBUG_MSG_REFRESH  
  CkPrintf ("DEBUG_CHARM %p update data_msg_ %p\n",this,data_msg_);
#endif
  if (data_msg_ == nullptr) return;

  data_msg_->update(data,is_local_);

#ifdef DEBUG_MSG_REFRESH  
  print(std::string(tag_)+":update");
#endif

  if (!is_local_) {
      CkFreeMsg (buffer_);
      buffer_ = nullptr;
  }
}

//----------------------------------------------------------------------

void MsgRefresh::print (std::string message)
{
#ifdef DEBUG_MSG_REFRESH  
  CkPrintf ("%s MSG_REFRESH %d %s %p %s %s\n",
            tag_,id_refresh_,message.c_str(),this,name_block_.c_str(),name_type_.c_str());
  if (data_msg_) {
    data_msg_->print((std::string(tag_) + ":" + message).c_str());
  } else {
    CkPrintf ("%s MSG_REFRESH data_msg_ = nil\n",tag_);
  }
#endif  
}
