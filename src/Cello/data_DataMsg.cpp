// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_DATA
// #define NEW_PARTICLE

//----------------------------------------------------------------------

int DataMsg::id_count = -1;

void * DataMsg::pack (DataMsg * msg)
{
#ifdef DEBUG_DATA
  cello::backtrace("pack");
#endif

  FieldFace    * ff = msg->field_face_;
  char         * fa = msg->field_array_;
#ifdef NEW_PARTICLE
  ParticleData * pd = msg->particle_data_;
  char         * pa = msg->particle_array_;
#endif

#ifdef DEBUG_DATA
  ff->print("packed");
#endif

#ifdef NEW_PARTICLE
  // Particle particle (particle_descr, msg->particle_data_);
#endif

  // Update ID (for debugging)

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();

  //--------------------------------------------------
  //  1. determine buffer size (must be consistent with #3)
  //--------------------------------------------------


#ifdef DEBUG_DATA
  CkPrintf ("%d %p DEBUG pack() %d\n",CkMyPe(),msg,msg->id_);
  fflush(stdout);
#endif

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
#ifdef NEW_PARTICLE
  ParticleDescr * particle_descr = simulation->particle_descr();
#endif

  Field field (field_descr, msg->field_data_);

  int size = 0;

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (fa) ? ff->num_bytes_array(field) : 0;
#ifdef NEW_PARTICLE
  const int n_pa = (pd) ? pd->data_size(particle_descr) : 0;
#ifdef DEBUG_DATA
  CkPrintf ("n_pa = %d\n",n_pa);
#endif
#endif

  size += sizeof(int); // n_ff_
  size += sizeof(int); // n_fa_
#ifdef NEW_PARTICLE
  size += sizeof(int); // n_pa_
#endif
  size += sizeof(int); // id_

  size += n_ff*sizeof(char);
  size += n_fa*sizeof(char);
#ifdef NEW_PARTICLE
  size += n_pa*sizeof(char);
#endif
  // (Particles)


  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

#ifdef DEBUG_DATA
  CkPrintf ("%d allocating buffer %p %d\n",__LINE__,msg,size); fflush(stdout);
#endif
  
  char * buffer = (char *) CkAllocBuffer (msg,size);
#ifdef DEBUG_DATA
  CkPrintf ("%d allocated buffer %p %d\n",__LINE__,msg,size); fflush(stdout);
#endif

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
#ifdef NEW_PARTICLE
  (*pi++) = n_pa;
#endif

  (*pi++) = msg->id_;

  if (n_ff > 0) {
    ff->save_data (pc);
    pc += n_ff;
  }
  if (n_fa > 0) {
    ff->face_to_array(field,pc);
    pc += n_fa;
  }
#ifdef NEW_PARTICLE
  if (n_pa > 0) {
    pd->save_data(particle_descr,pc);
    pc += n_pa;
  }
#endif

  // serialize ParticleData

  // ...

  //  CkFreeMsg(msg);
  delete msg;

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

DataMsg * DataMsg::unpack(void * buffer)
{
#ifdef DEBUG_DATA
  cello::backtrace("unpack");
#endif
  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  int size = sizeof(DataMsg);
#ifdef DEBUG_DATA
  CkPrintf ("%d allocating buffer %p %d\n",__LINE__,buffer,size); fflush(stdout);
#endif
  DataMsg * msg = (DataMsg *) CkAllocBuffer (buffer,sizeof(DataMsg));
  msg = new ((void*)msg) DataMsg;
#ifdef DEBUG_DATA
  CkPrintf ("%d allocated buffer %p %d\n",__LINE__,buffer,size); fflush(stdout);
#endif

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
#ifdef NEW_PARTICLE
  const int n_pa = (*pi++);
#endif
  msg->id_ = (*pi++);

#ifdef DEBUG_DATA
  CkPrintf ("%d %p DEBUG unpack() %d\n",CkMyPe(),msg,msg->id_);
  fflush(stdout);
#endif

  if (n_ff > 0) {
    msg->field_face_->load_data (pc);
    pc += n_ff;
  }

#ifdef DEBUG_DATA
  FieldFace * ff = msg->field_face_;
  ff->print("unpack");
#endif
  if (n_fa > 0) {
    msg->field_array_ = pc;
    pc += n_fa;
  } else {
    msg->field_array_ = NULL;
  }

#ifdef NEW_PARTICLE
  if (n_pa > 0) {
    msg->particle_array_ = pc;
    pc += n_pa;
  } else {
    msg->particle_array_ = NULL;
  }
#endif  

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data)
{
  
#ifdef DEBUG_DATA
  cello::backtrace("update");
  CkPrintf ("%d %p DEBUG update() %d\n",CkMyPe(),this,id_);
  fflush(stdout);
#endif

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  // Set field_face Field
  FieldDescr * field_descr = simulation->field_descr();
#ifdef NEW_PARTICLE
  ParticleDescr * particle_descr = simulation->particle_descr();
#endif
 
  Field field_dst = data->field();
 
  FieldFace * ff = field_face_;

  if (is_local_) {

#ifdef DEBUG_DATA
    ff->print("update local");
#endif
    //    field_face_->invert_face();

    Field field_src(field_descr,field_data_);
    field_face_->face_to_face(field_src, field_dst);

    delete field_face_;

  } else { // ! is_local_

#ifdef DEBUG_DATA    
    ff->print("update remote");
#endif

    if (field_array_ != NULL) {

      // Invert face since incoming not outgoing
#ifdef DEBUG_DATA
      CkPrintf ("%d %s:%d DEBUG\n",CkMyPe(),__FILE__,__LINE__);
       fflush(stdout);
#endif
      field_face_->invert_face();
#ifdef DEBUG_DATA
      CkPrintf ("%d %s:%d DEBUG\n",CkMyPe(),__FILE__,__LINE__);
       fflush(stdout);
#endif
      // @@@@ BUG
#ifdef DEBUG_DATA
      CkPrintf ("%d %s:%d DEBUG field_array_ = %p\n",CkMyPe(),__FILE__,__LINE__,
		field_array_); fflush(stdout);
      CkPrintf ("%d %s:%d DEBUG field_array_[0] = %d\n",CkMyPe(),__FILE__,__LINE__,
		field_array_[0]); fflush(stdout);
#endif
      field_face_->array_to_face(field_array_,field_dst);

    }

#ifdef NEW_PARTICLE
    if (particle_array_ != NULL) {

      // Invert face since incoming not outgoing
#ifdef DEBUG_DATA
      CkPrintf ("%d %s:%d DEBUG\n",CkMyPe(),__FILE__,__LINE__);
      fflush(stdout);
#endif
       //      field_face_->array_to_face(field_array_,field_dst);

    }
#endif
    CkFreeMsg (buffer_);
  }
}
