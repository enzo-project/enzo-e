// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"
#include "charm_simulation.hpp"

#define DEBUG_DATA

//----------------------------------------------------------------------

int DataMsg::id_count = -1;

void * DataMsg::pack (DataMsg * msg)
{

  FieldFace * ff = msg->field_face_;
  char      * fa = msg->field_array_;

  //--------------------------------------------------
  //  1. determine buffer size (must be consistent with #3)
  //--------------------------------------------------

  // FieldFace and FieldData serialized data

  FieldDescr * field_descr = proxy_simulation.ckLocalBranch()->field_descr();
  Field field (field_descr, msg->field_data_);

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();
  
  int size = 0;

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (fa) ? ff->num_bytes_array(field) : 0;

  size += sizeof(int); // n_ff_
  size += sizeof(int); // n_fa_
  size += sizeof(int); // id_
  size += n_ff*sizeof(char);
  size += n_fa*sizeof(char);
  // (Particles)


  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

  char * buffer = (char *) CkAllocBuffer (msg,size);

  //--------------------------------------------------
  //  3. serialize message data into buffer 
  //     (must be consistent with #1)
  //--------------------------------------------------

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  (*pi++) = n_ff;
  (*pi++) = n_fa;
  (*pi++) = msg->id_;

  if (n_ff > 0) {
    ff->save_data (pc);
    pc += n_ff;
  }
  if (n_fa > 0) {
    ff->face_to_array(field,pc);
    pc += n_fa;
  }

  // serialize ParticleData

  // ...

  delete msg;

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

DataMsg * DataMsg::unpack(void * buffer)
{
  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  DataMsg * msg = (DataMsg *) CkAllocBuffer (buffer,sizeof(DataMsg));

  msg = new ((void*)msg) DataMsg;

  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message
 
  union {
    char * pc;
    int  * pi;
  };

  pc = (char *)buffer;

  msg->field_face_ = new FieldFace;

  const int n_ff = (*pi++);
  const int n_fa = (*pi++);
  msg->id_ = (*pi++);

  if (n_ff > 0) {
    msg->field_face_->load_data (pc);
    pc += n_ff;
  }
  if (n_fa > 0) {
    msg->field_array_ = pc;
    pc += n_fa;
  } else {
    msg->field_array_ = NULL;
  }

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data)
{
  FieldDescr * field_descr = proxy_simulation.ckLocalBranch()->field_descr();
 
  // Invert face since incoming not outgoing
  Field field_dst = data->field();
 
  FieldFace * ff = field_face_;

  // Invert face since incoming not outgoing


  if (is_local_) {

    Field field_src(field_descr,field_data_);
    field_face_->face_to_face(field_src, field_dst);

    delete field_face_;

  } else { // ! is_local_

    if (field_array_ != NULL) {

      // Invert face since incoming not outgoing
      field_face_->invert_face();
      field_face_->array_to_face(field_array_,field_dst);

      CkFreeMsg (buffer_);
    }
  }
}
