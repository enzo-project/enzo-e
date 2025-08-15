// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgState.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2023-04-18
/// @brief    [\ref Charm] Declaration of the MsgState Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

//----------------------------------------------------------------------

long MsgState::counter[CONFIG_NODE_SIZE] = { };

//----------------------------------------------------------------------

MsgState::MsgState()
  : CMessage_MsgState(),
    time_(0.0),
    dt_(0.0),
    cycle_(0),
    stop_(0)
{
  ++counter[cello::index_static()];
}

//----------------------------------------------------------------------

MsgState::~MsgState()
{
  --counter[cello::index_static()];
}

//----------------------------------------------------------------------

void MsgState::update (Simulation * simulation)
{
  simulation->state()->init(cycle_,time_,dt_,stop_);
}

//----------------------------------------------------------------------

void * MsgState::pack (MsgState * msg)
{
  int size = 0;

  SIZE_SCALAR_TYPE(size,double,time_);
  SIZE_SCALAR_TYPE(size,double,dt_);
  SIZE_SCALAR_TYPE(size,int,cycle_);
  SIZE_SCALAR_TYPE(size,int,stop_);

  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

  char * buffer = (char *) CkAllocBuffer (msg,size);

  //--------------------------------------------------
  //  3. serialize message data into buffer 
  //--------------------------------------------------

  char * pc = buffer;

  SAVE_SCALAR_TYPE(pc,double,msg->time_);
  SAVE_SCALAR_TYPE(pc,double,msg->dt_);
  SAVE_SCALAR_TYPE(pc,int,msg->cycle_);
  SAVE_SCALAR_TYPE(pc,int,msg->stop_);

  delete msg;

  // Return the buffer

  ASSERT2("MsgState::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgState * MsgState::unpack(void * buffer)
{

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.

  MsgState * msg = 
    (MsgState *) CkAllocBuffer (buffer,sizeof(MsgState));

  msg = new ((void*)msg) MsgState;

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  char * pc = (char *) buffer;

  LOAD_SCALAR_TYPE(pc,double,msg->time_);
  LOAD_SCALAR_TYPE(pc,double,msg->dt_);
  LOAD_SCALAR_TYPE(pc,int,msg->cycle_);
  LOAD_SCALAR_TYPE(pc,int,msg->stop_);

  return msg;
}
