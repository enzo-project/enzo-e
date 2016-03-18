// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgCoarsen.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgCoarsen Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_DATA

//----------------------------------------------------------------------

int MsgCoarsen::id_count = -1;
long MsgCoarsen::counter = 0;

//----------------------------------------------------------------------

MsgCoarsen::MsgCoarsen()
    : CMessage_MsgCoarsen(),
      is_local_(true),
      id_(-1),
      data_msg_(new DataMsg),
      buffer_(NULL)
  {  ++counter; }

//----------------------------------------------------------------------

MsgCoarsen::~MsgCoarsen()
{
  --counter;
}

//----------------------------------------------------------------------

FieldFace * MsgCoarsen::field_face () 
{
  return (data_msg_) ? data_msg_->field_face() : NULL; 
}

//----------------------------------------------------------------------

void MsgCoarsen::set_data_msg  (DataMsg * data_msg) 
{
  if (data_msg_) {
    WARNING ("MsgCoarsen::set_data_msg()",
	     "overwriting existing data_msg_");
    delete data_msg_;
  }
  data_msg_ = data_msg;
}

//----------------------------------------------------------------------

// void MsgCoarsen::set_field_data  (FieldData * field_data) 
// {
//   data_msg_->set_field_data(field_data); 
// }

// //----------------------------------------------------------------------

// void MsgCoarsen::set_particle_data  (ParticleData * particle_data) 
// { 
//   data_msg_->set_particle_data(particle_data); 
// }

// //----------------------------------------------------------------------

// void MsgCoarsen::set_field_face    (FieldFace * field_face) 
// { 
//   data_msg_->set_field_face(field_face); 
// }

// //----------------------------------------------------------------------

void * MsgCoarsen::pack (MsgCoarsen * msg)
{

#ifdef DEBUG_DATA
  CkPrintf ("DEBUG %p MsgCoarsen::pack()\n",msg);
#endif
  // Update ID (for debugging)

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();

  int size = 0;

  size += sizeof(int); // id_
  size += msg->data_msg_->data_size();

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
  pc = msg->data_msg_->save_data(pc);

  delete msg;

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgCoarsen * MsgCoarsen::unpack(void * buffer)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  MsgCoarsen * msg = 
    (MsgCoarsen *) CkAllocBuffer (buffer,sizeof(MsgCoarsen));

  msg = new ((void*)msg) MsgCoarsen;
#ifdef DEBUG_DATA
  CkPrintf ("DEBUG %p MsgCoarsen::unpack()\n",msg);
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

  pc = msg->data_msg_->load_data(pc);

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void MsgCoarsen::update (Data * data)
{
#ifdef DEBUG_DATA
  CkPrintf ("DEBUG %p MsgCoarsen::update()\n",this);
#endif
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr    *    field_descr = simulation->   field_descr();
 
  Field field_dst = data->field();
 
  FieldData    * fd = data_msg_->field_data();
  ParticleData * pd = data_msg_->particle_data();
  FieldFace    * ff = data_msg_->field_face();
  char         * fa = data_msg_->field_array();

  if (pd != NULL) {

    // Insert new particles 

    Particle particle = data->particle();
    
    for (int it=0; it<particle.num_types(); it++) {
#ifdef DEBUG_DATA
      CkPrintf ("%d %p DEBUG update\n",CkMyPe(),pd);
      fflush(stdout);
#endif
      particle.gather (it, 1, &pd);
    }
    data_msg_->delete_particle_data();
  }

  if (fa != NULL) {

    if (is_local_) {

      Field field_src(field_descr,fd);
      ff->face_to_face(field_src, field_dst);

      delete ff;

    } else { // ! is_local_


      // Invert face since incoming not outgoing

      ff->invert_face();

      ff->array_to_face(fa,field_dst);


    }
    // @@@ BUG: segfaults
    //    delete field_face_;
    //    field_face_ = 0;
  }
  if (!is_local_) {
      CkFreeMsg (buffer_);
  }
}
