// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgCoarsen.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgCoarsen Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_NEW_REFRESH

//----------------------------------------------------------------------

int MsgCoarsen::id_count = -1;
long MsgCoarsen::counter = 0;

//----------------------------------------------------------------------

MsgCoarsen::MsgCoarsen()
  : CMessage_MsgCoarsen(),
    is_local_(true),
    id_(-1),
    data_msg_(NULL),
    buffer_(NULL),
    num_face_level_(0),
    face_level_(NULL)
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d MsgCoarsen::MsgCoarsen() %p\n",CkMyPe(),this);
  fflush(stdout);
#endif
  ic3_[0] = ic3_[1] = ic3_[2] = -1;
  ++counter; 
}

//----------------------------------------------------------------------

MsgCoarsen::MsgCoarsen(int num_face_level, int face_level[], int ic3[3])
  : CMessage_MsgCoarsen(),
    is_local_(true),
    id_(-1),
    data_msg_(NULL),
    buffer_(NULL),
    num_face_level_(num_face_level),
    face_level_(new int[num_face_level])
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d MsgCoarsen::MsgCoarsen() %p\n",CkMyPe(),this);
  fflush(stdout);
#endif

  ++counter; 

  for (int i=0; i<num_face_level_; i++) {
    face_level_[i] = face_level[i];
  }
  ic3_[0]=ic3[0];
  ic3_[1]=ic3[1];
  ic3_[2]=ic3[2];
}

//----------------------------------------------------------------------

MsgCoarsen::~MsgCoarsen()
{
  --counter;
  // delete data_msg_;
  // data_msg_ = 0;
  delete [] face_level_;
  face_level_ = 0;
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

void * MsgCoarsen::pack (MsgCoarsen * msg)
{

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgCoarsen::pack()\n",msg);
#endif
  // Update ID (for debugging)

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();

  int size = 0;

  // id_
  size += sizeof(int); 

  // have_data
  size += sizeof(int); 

  int have_data = (msg->data_msg_ != NULL);
  if (have_data) {
    size += msg->data_msg_->data_size();
  }

  // num_face_level_
  size += sizeof(int);

  // face_level_[]
  size += msg->num_face_level_ * sizeof(int);

  // ic3_[]
  size += 3*sizeof(int);

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

  // id_
  (*pi++) = msg->id_;

  have_data = (msg->data_msg_ != NULL);

  // have_data
  (*pi++) = have_data; 

  if (have_data) {
    // data_msg_
    pc = msg->data_msg_->save_data(pc);   
  }

  // num_face_level_
  (*pi++) = msg->num_face_level_;

  // face_level_[]
  for (size_t i=0; i<msg->num_face_level_; i++) {
    (*pi++) = msg->face_level_[i];
  }

  // ic3_[]
  (*pi++) = msg->ic3_[0];
  (*pi++) = msg->ic3_[1];
  (*pi++) = msg->ic3_[2];

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p MsgCoarsen pack message size %d %d\n",msg,(pc - (char*)buffer),size);
#endif

  ASSERT2("MsgRefresh::pack()",
	  "buffer size mismatch %d allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

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
#ifdef DEBUG_NEW_REFRESH
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

  // id_
  msg->id_ = (*pi++);

  // have_data
  int have_data = (*pi++);

  if (have_data) {
    // data_msg_
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = NULL;
  }

  // num_face_level_
  msg->num_face_level_ = (*pi++);

  // face_level_[]
  if (msg->num_face_level_ > 0) {
    msg->face_level_ = new int [msg->num_face_level_];
    for (size_t i = 0; i<msg->num_face_level_; i++) {
      msg->face_level_[i] = (*pi++);
    }
  } else {
    msg->face_level_ = 0;
  }

  // ic3_[]
  msg->ic3_[0] = (*pi++);
  msg->ic3_[1] = (*pi++);
  msg->ic3_[2] = (*pi++);

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void MsgCoarsen::update (Data * data)
{

  if (data_msg_ == NULL) return;

#ifdef DEBUG_NEW_REFRESH
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
#ifdef DEBUG_NEW_REFRESH
      CkPrintf ("%d %p DEBUG update\n",CkMyPe(),pd);
      fflush(stdout);
#endif
      particle.gather (it, 1, &pd);
    }
    // Don't delete particle data if local--done by child Block::data_
    // destructor()
    if (!is_local_) {
      data_msg_->delete_particle_data();
    }
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
