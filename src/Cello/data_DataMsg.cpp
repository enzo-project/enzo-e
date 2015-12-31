// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include "data.hpp"

//----------------------------------------------------------------------

void * DataMsg::pack (DataMsg * msg)
{
  //  1. determine buffer size

  int size = 0;

  // FieldFace
  size += 1;
  if (msg->ff_) size += msg->ff_->data_size();

  // FieldData
  size += 1;
  if (msg->fd_) size += msg->fd_->data_size();

  // ParticleData
  size += 1;
  if (msg->pd_) size += msg->pd_->data_size();

  //  2. allocate buffer using CkAllocBuffer()

  char * buffer = (char *) CkAllocBuffer (msg,size);

  //  3. serialize message data into buffer (along with control info
  //  required for de-serializing it)

  char * p = buffer;

  // FieldFace
  (*p++) = msg->ff_ ? 1 : 0;
  if (msg->ff_) p = msg->ff_->save_data (p);

  // FieldData
  (*p++) = msg->fd_ ? 1 : 0;
  if (msg->fd_) p = msg->fd_->save_data (p);

  // ParticleData
  (*p++) = msg->pd_ ? 1 : 0;
  if (msg->pd_) p = msg->pd_->save_data (p);

  //  4. free resources occupied by the message, including the message
  //  itself

  delete msg;

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

DataMsg * DataMsg::unpack(void * buffer)
{
  printf ("DataMsg::unpack()\n");
  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.

  DataMsg * msg = (DataMsg *) CkAllocBuffer (buffer,sizeof(DataMsg));

  // 2. De-serialize message data from input buffer into the allocated
  // message

  char * p = (char *) buffer;

  if (*p++) {
    msg->ff_ = new FieldFace;
    msg->ff_->load_data(p);
  }

  if (*p++) {
    msg->fd_ = new FieldData;
    msg->fd_->load_data(p);
  }

  if (*p++) {
    msg->pd_ = new ParticleData;
    msg->pd_->load_data(p);
  }

  // 3. Free the input buffer using CkFreeMsg.

  CkFreeMsg (buffer);

  return msg;
}

//----------------------------------------------------------------------
