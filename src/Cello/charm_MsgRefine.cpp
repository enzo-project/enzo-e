// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgRefine Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_MSG_REFINE

//----------------------------------------------------------------------

long MsgRefine::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

MsgRefine::MsgRefine()
  : CMessage_MsgRefine(),
    is_local_(true),
    data_msg_(NULL),
    buffer_(NULL),
    time_(-1.0), dt_(-1.0),
    index_(),
    nx_(-1), ny_(-1), nz_(-1),
    num_field_blocks_(-1),
    num_adapt_steps_(-1),
    cycle_(-1),
    refresh_type_(refresh_unknown),
    num_face_level_(0), face_level_(NULL)
{
  ++counter[cello::index_static()]; 
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE creating %p\n",CkMyPe(),__FILE__,__LINE__,this);
#endif  
}

//----------------------------------------------------------------------

MsgRefine::MsgRefine
(Index index,
 int nx, int ny, int nz,
 int num_field_blocks, int num_adapt_steps,
 int cycle, double time, double dt, int refresh_type,
 int num_face_level, int * face_level) 
  : CMessage_MsgRefine(),
    is_local_(true),
    data_msg_(NULL),
    buffer_(NULL),
    time_(time), dt_(dt),
    index_(index),
    nx_(nx), ny_(ny), nz_(nz),
    num_field_blocks_(num_field_blocks),
    num_adapt_steps_(num_adapt_steps),
    cycle_(cycle),
    refresh_type_(refresh_type),
  num_face_level_(num_face_level),
  face_level_(new int[num_face_level])
{  
  ++counter[cello::index_static()]; 
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE creating %p\n",CkMyPe(),__FILE__,__LINE__,this);
#endif  

  for (int i=0; i<num_face_level_; i++) {
    face_level_[i] = face_level[i];
  }
}

//----------------------------------------------------------------------

MsgRefine::~MsgRefine()
{
  --counter[cello::index_static()];
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE destroying %p\n",CkMyPe(),__FILE__,__LINE__,this);
#endif  

  delete data_msg_;
  data_msg_ = 0;
  delete [] face_level_;
  face_level_ = 0;
}

//----------------------------------------------------------------------

void MsgRefine::set_data_msg  (DataMsg * data_msg) 
{
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %p %s:%d DEBUG_MSG_REFINE MsgRefine::set_data_msg()\n",
	    CkMyPe(),this,__FILE__,__LINE__);
  fflush(stdout);
#endif  
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

#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE packing %p\n",CkMyPe(),__FILE__,__LINE__,msg);
#endif  

  // WARNING("MsgRefine::pack()",
  // 	  "message already has a buffer allocated");
  
  if (msg->buffer_ != NULL) return msg->buffer_;

  int size = 0;

  // ...determine buffer size

  size += sizeof(double);   // attribute-01
  size += sizeof(double);   // attribute-02
  size += 3*sizeof(int);   // attribute-03
  size += 3*sizeof(int);   // attribute-04
  size += sizeof(int);     // attribute-05
  size += sizeof(int);     // attribute-06
  size += sizeof(int);     // attribute-07
  size += sizeof(int);   // attribute-08
  size += sizeof(int);   // attribute-09
  size += msg->num_face_level_ * sizeof(int); // attribute-10

  size += sizeof(int);  // element-11

  int have_data = (msg->data_msg_ != NULL);
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d DEBUG_MSG_REFINE %s:%d msg->data_msg_ = %p\n",CkMyPe(),__FILE__,__LINE__,msg->data_msg_);
  fflush(stdout);
#endif  
  
  if (have_data) {
    size += msg->data_msg_->data_size(); // element-12
  }

  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

  char * buffer = (char *) CkAllocBuffer (msg,size);

#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE ENTER MsgRefine::pack() msg %p --> buffer %p\n",
	    CkMyPe(),__FILE__,__LINE__,msg,buffer);
  fflush(stdout);
#endif  
  
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

  (*pd++) = msg->time_;   // attribute-01
  (*pd++) = msg->dt_;     // attribute-02
  int v3[3];
  msg->index_.values(v3);
  (*pi++) = v3[0];        // attribute-03
  (*pi++) = v3[1];        // attribute-03
  (*pi++) = v3[2];        // attribute-03

  (*pi++) = msg->nx_;     // attribute-04  
  (*pi++) = msg->ny_;     // attribute-04
  (*pi++) = msg->nz_;     // attribute-04

  (*pi++) = msg->num_field_blocks_; // attribute-05
  (*pi++) = msg->num_adapt_steps_;  // attribute-06

  (*pi++) = msg->cycle_;  // attribute-07
  (*pi++) = msg->refresh_type_; // attribute-08
  (*pi++) = msg->num_face_level_; // attribute-09
  for (int i=0; i<msg->num_face_level_; i++) {
    (*pi++) = msg->face_level_[i]; // attribute-10
  }

  // data_msg_
  have_data = (msg->data_msg_ != NULL);
  (*pi++) = have_data;  // element-11

  if (have_data) {
    // data_msg_
    pc = msg->data_msg_->save_data(pc); // element-12
  }

  ASSERT2("MsgRefine::pack()",
	  "buffer size mismatch %d allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

#ifdef DEBUG_MSG_REFINE  
  //  CkPrintf ("%s:%d DEBUG_MSG_REFINE CkFreeMsg (%p)\n",__FILE__,__LINE__,msg);
#endif
  CkFreeMsg (msg);
  //  CkFreeMsg(msg);

  // Return the buffer

  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgRefine * MsgRefine::unpack(void * buffer)
{

  // 1. Allocate message using CkAllocBuffer.  NOTE do not use new.

  MsgRefine * msg =
    (MsgRefine *) CkAllocBuffer (buffer,sizeof(MsgRefine));

  msg = new ((void*)msg) MsgRefine;
  
  msg->is_local_ = false;

#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE unpacking %p\n",CkMyPe(),__FILE__,__LINE__,msg);
#endif  

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  union {
    double * pd;
    int    * pi;
    bool   * pb;
    char   * pc;
  };

  pc = (char *) buffer;

  msg->time_ = (*pd++);   // attribute-01
  msg->dt_   = (*pd++);   // attribute-02
  int v3[3];
  v3[0] = (*pi++);         // attribute-03
  v3[1] = (*pi++);         // attribute-03
  v3[2] = (*pi++);         // attribute-03
  msg->index_.set_values(v3);

  msg->nx_ = (*pi++);  // attribute-04
  msg->ny_ = (*pi++);  // attribute-04
  msg->nz_ = (*pi++);  // attribute-04

  msg->num_field_blocks_ = (*pi++); // attribute-05
  msg->num_adapt_steps_ = (*pi++);  // attribute-06
  msg->cycle_ = (*pi++);            // attribute-07
  msg->refresh_type_ = (*pi++);     // attribute-08
  msg->num_face_level_ = (*pi++);   // attribute-09

  if (msg->num_face_level_ > 0) {
    msg->face_level_ = new int [msg->num_face_level_];
    for (int i = 0; i<msg->num_face_level_; i++) {
      msg->face_level_[i] = (*pi++); // attribute-10
    }
  } else {
    msg->face_level_ = NULL;
 }
  
  int have_data = (*pi++);   // element-11

  if (have_data) {
    // data_msg_
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc); // element-12
  } else {
    msg->data_msg_ = NULL;
  }

  // 3. Save the input buffer for freeing later

#ifdef DEBUG_MSG_REFINE  
  //  CkPrintf ("%s:%d DEBUG_MSG_REFINE CkFreeMsg (%p)\n",__FILE__,__LINE__,buffer);
#endif  
  msg->buffer_ = buffer;
  //  CkFreeMsg(buffer);
  
  return msg;
}

//----------------------------------------------------------------------

void MsgRefine::update (Data * data)
{
  
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d %s:%d DEBUG_MSG_REFINE updating %p\n",CkMyPe(),__FILE__,__LINE__,this);
#endif  

  if (data_msg_ == NULL) return;

  Simulation * simulation  = cello::simulation();
  FieldDescr * field_descr = cello::field_descr();
 
  Field field_dst = data->field();
 
  FieldData    * fd = data_msg_->field_data();
  ParticleData * pd = data_msg_->particle_data();
  FieldFace    * ff = data_msg_->field_face();
  char         * fa = data_msg_->field_array();

  if (pd != NULL) {

    // Insert new particles 

    Particle particle = data->particle();

    int count = 0;
    for (int it=0; it<particle.num_types(); it++) {
      count += particle.gather (it, 1, &pd);
    }
    simulation->data_insert_particles(count);

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
  }
  if (! is_local_) CkFreeMsg (buffer_);
}
