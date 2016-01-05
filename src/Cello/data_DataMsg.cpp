// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include "data.hpp"

//----------------------------------------------------------------------

void * DataMsg::pack (DataMsg * msg)
{

  CkPrintf ("DataMsg::unpack()\n");

  //--------------------------------------------------
  //  1. determine buffer size
  //--------------------------------------------------

  int size = 0;

  // FieldFace num_bytes
  size += sizeof(int);

  // FieldFace
  if (msg->field_face_) size += msg->field_face_->data_size();

  // FieldData
  size += sizeof(int);

  int nf = msg->field_data_ ? msg->field_face_->num_bytes_array() : 0;

  size += nf;

  // // ParticleData
  // size += 1;
  // if (msg->pd_) size += msg->pd_->data_size();

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

  // FieldFace
  (*pc++) = msg->field_face_ ? 1 : 0;

  if (msg->field_face_) {
    pc = msg->field_face_->save_data (pc);

  // FieldData buffer length

    int nf = msg->field_face_->num_bytes_array();

    (*pi++) = nf;
  
  // FieldData buffer data

    if (nf > 0) {
      msg->field_face_->face_to_array(pc);
      pc += nf;
    }
  }

  // // ParticleData
  // 
  // if (msg->pd_) p = msg->pd_->save_data (p);

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

  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message
 
  // char * p = (char *) buffer;
  // if (*p++) {
  //   msg->ff_ = new FieldFace;
  //   msg->ff_->load_data(p);
  //  }

  // if (*p++) {
  //   msg->fd_ = new FieldData;
  //   msg->fd_->load_data(p);
  // }

  // if (*p++) {
  //   msg->pd_ = new ParticleData;
  //   msg->pd_->load_data(p);
  // }

  // // 3. Free the input buffer using CkFreeMsg.

  CkFreeMsg (buffer);

  return msg;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data)
{
  if (is_local_) {

    // copy portion of field_data_ described by FieldFace to data
    union {
      char * pc;
      int  * pi;
    };

    pc = field_array_;

    // FieldFace
    (*pc++) = field_face_ ? 1 : 0;
    if (field_face_) pc = field_face_->save_data (pc);

    // FieldData buffer length

    int nf = field_face_->num_bytes_array();
    (*pi++) = nf;
  
    // FieldData buffer data

    if (nf > 0) {
      field_face_->face_to_array(pc);
      pc += nf;
    }


  } else {
    // copy field_array_ to data field
    
  }
}
