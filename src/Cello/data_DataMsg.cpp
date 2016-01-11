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

  msg->n_ff_ = (ff) ? ff->data_size() : 0;
  msg->n_fa_ = (fa) ? ff->num_bytes_array() : 0;

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();
  
  int size = 0;

  size += sizeof(int); // n_ff_
  size += sizeof(int); // n_fa_
  size += sizeof(int); // id_
  size += msg->n_ff_*sizeof(char);
  size += msg->n_fa_*sizeof(char);
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

  (*pi++) = msg->n_ff_;
  (*pi++) = msg->n_fa_;
  (*pi++) = msg->id_;

  if (msg->n_ff_ > 0) {
    ff->save_data (pc);
    pc += msg->n_ff_;
  }
  if (msg->n_fa_ > 0) {
    ff->face_to_array(pc);
    pc += msg->n_fa_;
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

  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message
 
  union {
    char * pc;
    int  * pi;
  };

  pc = (char *)buffer;

  msg->field_face_ = new FieldFace;

  msg->n_ff_ = (*pi++);
  msg->n_fa_ = (*pi++);
  msg->id_ = (*pi++);

  if (msg->n_ff_ > 0) {
    msg->field_face_->load_data (pc);
    pc += msg->n_ff_;
  }
  if (msg->n_fa_ > 0) {
    msg->field_array_ = pc;
    pc += msg->n_fa_;
  }

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data)
{
  
  int n = 0;
  char * array = NULL;

  if (is_local_) {
    field_face_->face_to_array(&n, &array);
  }

  // Set field_face Field
  Field field = field_face_->field();
  field.set_field_descr(proxy_simulation.ckLocalBranch()->field_descr());
  field.set_field_data (data->field_data());
  field_face_->set_field(field);

  // Invert face since incoming not outgoing

  field_face_->invert_face();

  if (is_local_) {

    field_face_->array_to_face(n,array);

    double * darray = (double*)array;

    delete field_face_;
    delete [] array;
    //    field_face_->loop_limits 
  } else {
      double * darray = (double*)field_array_;
      field_face_->array_to_face(n_fa_,field_array_);


      CkFreeMsg (buffer_);
  }
}
