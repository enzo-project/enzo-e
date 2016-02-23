// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_DATA

//----------------------------------------------------------------------

int DataMsg::id_count = -1;

//----------------------------------------------------------------------

void * DataMsg::pack (DataMsg * msg)
{

  FieldFace    * ff = msg->field_face_;
  char         * fa = msg->field_array_;
  ParticleData * pd = msg->particle_data_;

  // Update ID (for debugging)

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();

  //--------------------------------------------------
  //  1. determine buffer size (must be consistent with #3)
  //--------------------------------------------------

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  ParticleDescr * particle_descr = simulation->particle_descr();

  Field field (field_descr, msg->field_data_);

  int size = 0;

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (fa) ? ff->num_bytes_array(field) : 0;
  const int n_pa = (pd) ? pd->data_size(particle_descr) : 0;

  size += sizeof(int); // n_ff_
  size += sizeof(int); // n_fa_
  size += sizeof(int); // n_pa_
  size += sizeof(int); // id_

  size += n_ff*sizeof(char);
  size += n_fa*sizeof(char);
  size += n_pa*sizeof(char);

  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

  char * buffer = (char *) CkAllocBuffer (msg,size);

  //--------------------------------------------------
  //  3. serialize message data into buffer 
  //     (must be consistent with #1 and unpack())
  //--------------------------------------------------

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  (*pi++) = n_ff;
  (*pi++) = n_fa;
  (*pi++) = n_pa;
  (*pi++) = msg->id_;

  if (n_ff > 0) {
    ff->save_data (pc);
    pc += n_ff;
  }
  if (n_fa > 0) {
    ff->face_to_array(field,pc);
    pc += n_fa;
  }
  if (n_pa > 0) {
    pd->save_data(particle_descr,pc);
    pc += n_pa;
  }

  delete msg;

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

DataMsg * DataMsg::unpack(void * buffer)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  ParticleDescr * particle_descr = simulation->particle_descr();

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  DataMsg * msg = (DataMsg *) CkAllocBuffer (buffer,sizeof(DataMsg));
  msg = new ((void*)msg) DataMsg;

  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())
 
  union {
    char * pc;
    int  * pi;
  };

  pc = (char *) buffer;

  msg->field_face_ = new FieldFace;

  const int n_ff = (*pi++);
  const int n_fa = (*pi++);
  const int n_pa = (*pi++);
  msg->id_ = (*pi++);

  if (n_ff > 0) {
    pc = msg->field_face_->load_data (pc);
  }
  if (n_fa > 0) {
    msg->field_array_ = pc;
    pc += n_fa;
  } else {
    msg->field_array_ = NULL;
  }
  if (n_pa > 0) {
    ParticleData * pd = msg->particle_data_ = new ParticleData;
#ifdef DEBUG_DATA
    CkPrintf ("%d %p DEBUG unpack\n",CkMyPe(),pd);
    fflush(stdout);
#endif
    pd->allocate(particle_descr);
    pc = pd->load_data(particle_descr,pc);
  } else {
    msg->particle_data_ = NULL;
  }

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr    *    field_descr = simulation->   field_descr();
  ParticleDescr * particle_descr = simulation->particle_descr();
 
  Field field_dst = data->field();
 
  FieldFace * ff = field_face_;

  if (particle_data_ != NULL) {

    // Insert new particles 

    Particle particle = data->particle();
    
    for (int it=0; it<particle.num_types(); it++) {
#ifdef DEBUG_DATA
      CkPrintf ("%d %p DEBUG update\n",CkMyPe(),particle_data_);
      fflush(stdout);
#endif
      particle.gather (it, 1, &particle_data_);
    }
    delete particle_data_;
    particle_data_ = NULL;
  }

  if (field_array_ != NULL) {

    if (is_local_) {

      Field field_src(field_descr,field_data_);
      ff->face_to_face(field_src, field_dst);

      delete ff;

    } else { // ! is_local_


      // Invert face since incoming not outgoing

      ff->invert_face();

      ff->array_to_face(field_array_,field_dst);


    }
  }
  if (!is_local_) {
      CkFreeMsg (buffer_);
  }
}
