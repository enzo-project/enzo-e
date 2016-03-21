// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgRefine Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_NEW_REFRESH

//----------------------------------------------------------------------

int MsgRefine::id_count = -1;
long MsgRefine::counter = 0;

//----------------------------------------------------------------------

MsgRefine::MsgRefine()
  : CMessage_MsgRefine(),
    is_local_(true),
    id_(-1),
    data_msg_(NULL),
    index_(),
    nx_(-1), ny_(-1), nz_(-1),
    num_field_blocks_(-1),
    num_adapt_steps_(-1),
    cycle_(-1), time_(-1.0), dt_(-1.0),
    refresh_type_(refresh_unknown),
    num_face_level_(0), face_level_(NULL),
    testing_(false) ,
    buffer_(NULL)
{  
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d MsgRefine::MsgRefine() %p\n",CkMyPe(),this);
  fflush(stdout);
#endif
  ++counter; 
}

//----------------------------------------------------------------------

MsgRefine::MsgRefine
(Index index,
 int nx, int ny, int nz,
 int num_field_blocks, int num_adapt_steps,
 int cycle, double time, double dt, int refresh_type,
 int num_face_level, int * face_level, bool testing ) 
  : CMessage_MsgRefine(),
    is_local_(true),
    id_(-1),
    data_msg_(NULL),
    index_(index),
    nx_(nx), ny_(ny), nz_(nz),
    num_field_blocks_(num_field_blocks),
    num_adapt_steps_(num_adapt_steps),
    cycle_(cycle), time_(time), dt_(dt),
    refresh_type_(refresh_type),
    num_face_level_(num_face_level), face_level_(new int[27]),
  testing_(testing) ,
  buffer_(NULL)
{  
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d MsgRefine::MsgRefine() %p\n",CkMyPe(),this);
  CkPrintf ("%d %d %d\n",nx,ny,nz);
  CkPrintf ("%d %d\n",num_field_blocks,num_adapt_steps);
  CkPrintf ("%d %f %f %d\n",cycle,time,dt,refresh_type);
  CkPrintf ("%d %p %d\n",num_face_level,face_level,testing);
  fflush(stdout);
#endif

  ++counter; 
  for (int i=0; i<num_face_level_; i++) {
    face_level_[i] = face_level[i];
  }
}

//----------------------------------------------------------------------

MsgRefine::~MsgRefine()
{
  --counter;
  delete [] face_level_;
  face_level_ = 0;
}

//----------------------------------------------------------------------

void MsgRefine::set_data_msg  (DataMsg * data_msg) 
{
  if (data_msg_) {
    WARNING ("MsgRefine::set_data_msg()",
	     "overwriting existing data_msg_");
    delete data_msg_;
  }
  data_msg_ = data_msg;
}

//----------------------------------------------------------------------

void * MsgRefine::pack (MsgRefine * msg)
{

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d DEBUG MsgRefine::pack()\n",CkMyPe());
#endif
  // Update ID (for debugging)

  if (id_count == -1) id_count = CkMyPe()+CkNumPes();

  msg->id_ = id_count;

  id_count += CkNumPes();

  int size = 0;

  // ...determine buffer size

  // time_;
  size += sizeof(double);

  // dt_;
  size += sizeof(double);

  // id_
  size += sizeof(int);  

  // index_
  size += 3*sizeof(int);

  // nx_, ny_, nz_
  size += 3*sizeof(int);

  // num_field_blocks_
  size += sizeof(int);  

  // num_adapt_steps_
  size += sizeof(int);  

  // cycle_
  size += sizeof(int);  

  // refresh_type_
  size += sizeof(int);

  // num_face_level_
  size += sizeof(int);

  // face_level_[]
  size += msg->num_face_level_ * sizeof(int);

  // testing_
  size += sizeof(bool);

  // data_msg_

  int have_data = (msg->data_msg_ != NULL);
  if (have_data) {
    size += msg->data_msg_->data_size();
  }

  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

  char * buffer = (char *) CkAllocBuffer (msg,size);

  //--------------------------------------------------
  //  3. serialize message data into buffer 
  //--------------------------------------------------

  union {
    double * pd;
    int    * pi;
    bool   * pb;
    char   * pc;
  };

  pc = buffer;

  // time_
  (*pd++) = msg->time_;

  // dt_
  (*pd++) = msg->dt_;

  // id_
  (*pi++) = msg->id_;

  // index_
  int v3[3];
  msg->index_.values(v3);
  (*pi++) = v3[0];
  (*pi++) = v3[1];
  (*pi++) = v3[2];

  // nx_, ny_, nz_
  (*pi++) = msg->nx_;
  (*pi++) = msg->ny_;
  (*pi++) = msg->nz_;

  // num_field_blocks_
  (*pi++) = msg->num_field_blocks_;

  // num_adapt_steps_
  (*pi++) = msg->num_adapt_steps_;

  // cycle_
  (*pi++) = msg->cycle_;

  // refresh_type_
  (*pi++) = msg->refresh_type_;

  // num_face_level_
  (*pi++) = msg->num_face_level_;

  // face_level_[]
  for (size_t i=0; i<msg->num_face_level_; i++) {
    (*pi++) = msg->face_level_[i];
  }

  // testing_
  (*pb++) = msg->testing_;

  // data_msg_
  have_data = (msg->data_msg_ != NULL);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p MsgRefine pack message size %d\n",msg,(pc - (char*)buffer));
#endif

  delete msg;

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgRefine * MsgRefine::unpack(void * buffer)
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d MsgRefine::unpack()\n",CkMyPe());
  fflush(stdout);
#endif

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.
 
  MsgRefine * msg = 
    (MsgRefine *) CkAllocBuffer (buffer,sizeof(MsgRefine));

  msg = new ((void*)msg) MsgRefine;
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgRefine::unpack()\n",msg);
#endif
  
  msg->is_local_ = false;

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  union {
    double * pd;
    int    * pi;
    bool   * pb;
    char   * pc;
  };

  pc = (char *) buffer;

  // time_;
  msg->time_ = (*pd++);

  // dt_;
  msg->dt_ = (*pd++);

  // id_
  msg->id_ = (*pi++);

  // index_
  int v3[3];
  v3[0] = (*pi++);
  v3[1] = (*pi++);
  v3[2] = (*pi++);
  msg->index_.set_values(v3);

  // nx_, ny_, nz_
  msg->nx_ = (*pi++);
  msg->ny_ = (*pi++);
  msg->nz_ = (*pi++);

  // num_field_blocks_
  msg->num_field_blocks_ = (*pi++);
  
  // num_adapt_steps_
  msg->num_adapt_steps_ = (*pi++);

  // cycle_
  msg->cycle_ = (*pi++);

  // refresh_type_
  msg->refresh_type_ = (*pi++);

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
  
  // testing_
  msg->testing_ = (*pb++);

  // data_msg_
  int have_data = (*pi++);
  if (have_data) {
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = 0;
  }

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p MsgRefine unpack message size %d\n",msg,(pc - (char*)buffer));
#endif

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void MsgRefine::update (Data * data)
{
  if (data_msg_ == NULL) return;

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("DEBUG %p MsgRefine::update()\n",this);
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
