// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include "data.hpp"

//----------------------------------------------------------------------

void * DataMsg::pack (DataMsg * dm)
{
  CkPrintf ("DataMsg::pack()\n");

  //  1. determine buffer size

  int size = 0;

  //  if (pd_) size += pd_->data_size();
  //  if (fd_) size += fd_->data_size();
  if (dm->ff_) size += dm->ff_->data_size();

  //  2. allocate buffer using CkAllocBuffer()
  CkAllocBuffer (dm,size);

  //  3. serialize message data into buffer (along with control info
  //  required for de-serializing it)

  //  4. free resources occupied by the message, including the message
  //  itself

  delete dm;

}

//----------------------------------------------------------------------

DataMsg * DataMsg::unpack(void * buffer)
{
  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.

  // 2. De-serialize message data from input buffer into the allocated
  // message

  // 3. Free the input buffer using CkFreeMsg.

  printf ("DataMsg::unpack()\n");
  DataMsg * msg = new DataMsg;
  CkFreeMsg (buffer);
  return msg;
}

//----------------------------------------------------------------------
