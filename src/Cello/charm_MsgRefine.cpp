// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgRefine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    [\ref Charm] Declaration of the MsgRefine Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

//----------------------------------------------------------------------

long MsgRefine::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

MsgRefine::MsgRefine()
  : CMessage_MsgRefine(),
    is_local_(true),
    data_msg_(nullptr),
    index_(),
    nx_(-1), ny_(-1), nz_(-1),
    num_field_blocks_(-1),
    num_adapt_steps_(-1),
    refresh_type_(refresh_unknown),
    face_level_(),
    adapt_parent_(nullptr),
    state_(nullptr),
    restart_io_reader_(-1),
    buffer_(nullptr)
{
  ++counter[cello::index_static()];
}

//----------------------------------------------------------------------

MsgRefine::MsgRefine
(Index index,
 int nx, int ny, int nz,
 int num_field_blocks, int num_adapt_steps,
 int refresh_type,
 const std::vector<int> & face_level,
 Adapt * adapt_parent,
 State * state,
 int restart_io_reader
 )
  : CMessage_MsgRefine(),
    is_local_(true),
    data_msg_(nullptr),
    index_(index),
    nx_(nx), ny_(ny), nz_(nz),
    num_field_blocks_(num_field_blocks),
    num_adapt_steps_(num_adapt_steps),
    refresh_type_(refresh_type),
    face_level_(),
    adapt_parent_(adapt_parent),
    state_(state),
    restart_io_reader_(restart_io_reader),
    buffer_(nullptr)
{
  ++counter[cello::index_static()];
  face_level_ = face_level;
}

//----------------------------------------------------------------------

MsgRefine::~MsgRefine()
{
  --counter[cello::index_static()];

  delete data_msg_;
  data_msg_ = 0;
  CkFreeMsg (buffer_);
  buffer_=nullptr;
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

  // WARNING("MsgRefine::pack()",
  // 	  "message already has a buffer allocated");

  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = 0;

  // ...determine buffer size

  SIZE_OBJECT_PTR_TYPE(size,DataMsg,msg->data_msg_);
  SIZE_OBJECT_TYPE(    size,        msg->index_);
  SIZE_SCALAR_TYPE(    size,int,    msg->nx_);
  SIZE_SCALAR_TYPE(    size,int,    msg->ny_);
  SIZE_SCALAR_TYPE(    size,int,    msg->nz_);
  SIZE_SCALAR_TYPE(    size,int,    msg->num_field_blocks_);
  SIZE_SCALAR_TYPE(    size,int,    msg->num_adapt_steps_);
  SIZE_SCALAR_TYPE(    size,int,    msg->refresh_type_);
  SIZE_VECTOR_TYPE(    size,int,    msg->face_level_);
  SIZE_OBJECT_PTR_TYPE(size,Adapt,  msg->adapt_parent_);
  SIZE_OBJECT_PTR_TYPE(size,State,  msg->state_);
  SIZE_SCALAR_TYPE(    size,int,    msg->restart_io_reader_);

  //--------------------------------------------------
  //  2. allocate buffer using CkAllocBuffer()
  //--------------------------------------------------

  char * buffer = (char *) CkAllocBuffer (msg,size);

  //--------------------------------------------------
  //  3. serialize message data into buffer
  //--------------------------------------------------

  char * pc = (char *) buffer;

  SAVE_OBJECT_PTR_TYPE(pc,DataMsg,msg->data_msg_);
  SAVE_OBJECT_TYPE(    pc,        msg->index_);
  SAVE_SCALAR_TYPE(    pc,int,    msg->nx_);
  SAVE_SCALAR_TYPE(    pc,int,    msg->ny_);
  SAVE_SCALAR_TYPE(    pc,int,    msg->nz_);
  SAVE_SCALAR_TYPE(    pc,int,    msg->num_field_blocks_);
  SAVE_SCALAR_TYPE(    pc,int,    msg->num_adapt_steps_);
  SAVE_SCALAR_TYPE(    pc,int,    msg->refresh_type_);
  SAVE_VECTOR_TYPE(    pc,int,    msg->face_level_);
  SAVE_OBJECT_PTR_TYPE(pc,Adapt,  msg->adapt_parent_);
  SAVE_OBJECT_PTR_TYPE(pc,State,  msg->state_);
  SAVE_SCALAR_TYPE(    pc,int,    msg->restart_io_reader_);

  ASSERT2("MsgRefine::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  delete msg;

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

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())

  char * pc = (char *) buffer;

  LOAD_OBJECT_PTR_TYPE(pc,DataMsg,msg->data_msg_);
  LOAD_OBJECT_TYPE(    pc,        msg->index_);
  LOAD_SCALAR_TYPE(    pc,int,    msg->nx_);
  LOAD_SCALAR_TYPE(    pc,int,    msg->ny_);
  LOAD_SCALAR_TYPE(    pc,int,    msg->nz_);
  LOAD_SCALAR_TYPE(    pc,int,    msg->num_field_blocks_);
  LOAD_SCALAR_TYPE(    pc,int,    msg->num_adapt_steps_);
  LOAD_SCALAR_TYPE(    pc,int,    msg->refresh_type_);
  LOAD_VECTOR_TYPE(    pc,int,    msg->face_level_);
  LOAD_OBJECT_PTR_TYPE(pc,Adapt,  msg->adapt_parent_);
  LOAD_OBJECT_PTR_TYPE(pc,State,  msg->state_);
  LOAD_SCALAR_TYPE(    pc,int,    msg->restart_io_reader_);

  // 3. Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void MsgRefine::update (Data * data)
{

  if (data_msg_ == nullptr) return;

  Simulation * simulation  = cello::simulation();
  FieldDescr * field_descr = cello::field_descr();

  Field field_dst = data->field();

  FieldData    * fd = data_msg_->field_data();
  ParticleData * pd = data_msg_->particle_data();
  FieldFace    * ff = data_msg_->field_face();
  char         * fa = data_msg_->field_array();

  if (pd != nullptr) {

    // Insert new particles

    Particle particle = data->particle();

    int count = 0;
    for (int it=0; it<particle.num_types(); it++) {
      count += particle.gather (it, 1, &pd);
    }
    simulation->data_insert_particles(count);

    data_msg_->delete_particle_data();
  }

  if (fa != nullptr) {

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

  // Update scalar data
  data_msg_->update_scalars(data);

  if (! is_local_) {
    CkFreeMsg (buffer_);
    buffer_ = nullptr;
  }
}

//----------------------------------------------------------------------

void MsgRefine::print()
{
  CkPrintf ("bool is_local_ = %d\n", is_local_);
  CkPrintf ("DataMsg * data_msg_ = %p\n", (void*)data_msg_);
  CkPrintf ("void * buffer_ = %p\n",buffer_);
  CkPrintf ("\n");
  CkPrintf ("int nx_, ny_, nz_ = %d %d %d\n",nx_, ny_, nz_);
  CkPrintf ("int num_field_blocks_ = %d\n",num_field_blocks_);
  CkPrintf ("int num_adapt_steps_ = %d\n",num_adapt_steps_);
  CkPrintf ("int refresh_type_ = %d\n",refresh_type_);
  if (restart_io_reader_ >= 0)
    CkPrintf ("int restart_io_reader_ = %d\n",restart_io_reader_);
  if (adapt_parent_)
    adapt_parent_->print("MsgRefine");
  else
    CkPrintf ("Adapt * adapt_parent_ = nullptr\n");
}
