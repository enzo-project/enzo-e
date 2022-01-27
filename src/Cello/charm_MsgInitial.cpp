// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgInitial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-05-31
/// @brief    [\ref Charm] Declaration of the MsgInitial Charm++ message

#include "data.hpp"
#include "charm.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_MSG_INITIAL

#ifdef DEBUG_MSG_INITIAL
#  define TRACE_MSG_INITIAL(MSG) CkPrintf ("TRACE_MSG_INITIAL %d %s %d\n", \
                                           CkMyPe(),MSG,__LINE__); fflush(stdout);
#else
#  define TRACE_MSG_INITIAL(MSG) /* ... */
#endif

//----------------------------------------------------------------------

long MsgInitial::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

MsgInitial::MsgInitial()
  : CMessage_MsgInitial(),
    is_local_(true),
    buffer_(nullptr),
    data_type_(),
    data_name_(),
    data_attribute_(),
    data_precision_(),
    data_bytes_(0),
    data_values_(nullptr),
    data_delete_(true),
    count_(false),
    tag_(),
    n4_(),
    h4_(),
    nx_(0),ny_(0),nz_(0),
    IX_(0),IY_(0),IZ_(0)
{
  ++counter[cello::index_static()];
  cello::hex_string(tag_,TAG_LEN);
}

//----------------------------------------------------------------------

MsgInitial::~MsgInitial()
{
  --counter[cello::index_static()];
  if (data_delete_) delete [] data_values_;
  CkFreeMsg (buffer_);
  buffer_=nullptr;
}

//----------------------------------------------------------------------

void * MsgInitial::pack (MsgInitial * msg)
{
  TRACE_MSG_INITIAL("pack");
  // Return with buffer if already packed
  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = 0;

  // determine buffer size

  SIZE_STRING_TYPE(size,msg->data_type_);
  SIZE_STRING_TYPE(size,msg->data_name_);
  SIZE_STRING_TYPE(size,msg->data_attribute_);
  SIZE_SCALAR_TYPE(size,int, msg->data_precision_);
  SIZE_SCALAR_TYPE(size,int, msg->data_bytes_);
  SIZE_ARRAY_TYPE (size,char,msg->data_values_,msg->data_bytes_);
  SIZE_SCALAR_TYPE(size,int, msg->data_delete_);
  SIZE_SCALAR_TYPE(size,int, msg->count_);
  SIZE_ARRAY_TYPE (size,char,msg->tag_,TAG_LEN+1);
  SIZE_ARRAY_TYPE (size,int,   msg->n4_,4);
  SIZE_ARRAY_TYPE (size,double,msg->h4_,4);
  SIZE_SCALAR_TYPE(size,int, msg->nx_);
  SIZE_SCALAR_TYPE(size,int, msg->ny_);
  SIZE_SCALAR_TYPE(size,int, msg->nz_);
  SIZE_SCALAR_TYPE(size,int, msg->IX_);
  SIZE_SCALAR_TYPE(size,int, msg->IY_);
  SIZE_SCALAR_TYPE(size,int, msg->IZ_);

  //--------------------------------------------------

  // allocate buffer using CkAllocBuffer()

  char * buffer = (char *) CkAllocBuffer (msg,size);

  // serialize message data into buffer 

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  SAVE_STRING_TYPE(pc,msg->data_type_);
  SAVE_STRING_TYPE(pc,msg->data_name_);
  SAVE_STRING_TYPE(pc,msg->data_attribute_);
  SAVE_SCALAR_TYPE(pc,int,msg->data_precision_);
  SAVE_SCALAR_TYPE(pc,int,msg->data_bytes_);
  SAVE_ARRAY_TYPE(pc,char,msg->data_values_,msg->data_bytes_);
  SAVE_SCALAR_TYPE(pc,int,msg->data_delete_);
  SAVE_SCALAR_TYPE(pc,int,msg->count_);
  SAVE_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);
  SAVE_ARRAY_TYPE (pc,int,   msg->n4_,4);
  SAVE_ARRAY_TYPE (pc,double,msg->h4_,4);
  SAVE_SCALAR_TYPE(pc,int, msg->nx_);
  SAVE_SCALAR_TYPE(pc,int, msg->ny_);
  SAVE_SCALAR_TYPE(pc,int, msg->nz_);
  SAVE_SCALAR_TYPE(pc,int, msg->IX_);
  SAVE_SCALAR_TYPE(pc,int, msg->IY_);
  SAVE_SCALAR_TYPE(pc,int, msg->IZ_);

  ASSERT2("MsgInitial::pack()",
          "buffer size mismatch %ld allocated %d packed",
          (pc - (char*)buffer),size,
          (pc - (char*)buffer) == size);

  CkFreeMsg (msg);
  // Return the buffer
  return (void *) buffer;
}

//----------------------------------------------------------------------

MsgInitial * MsgInitial::unpack(void * buffer)
{

  TRACE_MSG_INITIAL("unpack");
  // Allocate storage using CkAllocBuffer (not new!)

  MsgInitial * msg = (MsgInitial *) CkAllocBuffer (buffer,sizeof(MsgInitial));

  msg = new ((void*)msg) MsgInitial;

  msg->is_local_ = false;

  // de-serialize message data from input buffer into allocated message

  union {
    char   * pc;
    int    * pi;
  };

  pc = (char *) buffer;

  LOAD_STRING_TYPE(pc,msg->data_type_);
  LOAD_STRING_TYPE(pc,msg->data_name_);
  LOAD_STRING_TYPE(pc,msg->data_attribute_);
  LOAD_SCALAR_TYPE(pc,int,msg->data_precision_);
  LOAD_SCALAR_TYPE(pc,int,msg->data_bytes_);
  msg->data_values_ = new char[msg->data_bytes_];
  LOAD_ARRAY_TYPE (pc,char,msg->data_values_,msg->data_bytes_);
  LOAD_SCALAR_TYPE(pc,int,msg->data_delete_);
  msg->data_delete_ = true;
  LOAD_SCALAR_TYPE(pc,int,msg->count_);
  LOAD_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);
  LOAD_ARRAY_TYPE (pc,int,   msg->n4_,4);
  LOAD_ARRAY_TYPE (pc,double,msg->h4_,4);
  LOAD_SCALAR_TYPE(pc,int, msg->nx_);
  LOAD_SCALAR_TYPE(pc,int, msg->ny_);
  LOAD_SCALAR_TYPE(pc,int, msg->nz_);
  LOAD_SCALAR_TYPE(pc,int, msg->IX_);
  LOAD_SCALAR_TYPE(pc,int, msg->IY_);
  LOAD_SCALAR_TYPE(pc,int, msg->IZ_);

  // Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void MsgInitial::update (Data * data)
{
  TRACE_MSG_INITIAL("update");
  // Copy field or particle data into data
  if (!is_local_) {
    CkFreeMsg (buffer_);
    buffer_ = nullptr;
  } 
}

//----------------------------------------------------------------------

void MsgInitial::print (const char * msg)
{
  CkPrintf ("MSG_INITIAL====================\n");
  CkPrintf ("MSG_INITIAL tag %s %s\n",msg,tag_);
  CkPrintf ("MSG_INITIAL is_local %d\n",is_local_);
}
//----------------------------------------------------------------------

void MsgInitial::set_field_data
(std::string field_name,
 char * data, int data_size, int data_precision)
{
  data_type_      = "field";
  data_name_      = field_name;
  data_attribute_ = "";  // unused for field data
  
  copy_data_(data,data_size,data_precision);
}


//----------------------------------------------------------------------

void MsgInitial::get_field_data
(std::string * field_name,
 char ** data, int * data_precision)
{
  (*field_name)     = data_name_;
  (*data)           = data_values_;
  (*data_precision) = data_precision_;
}


//----------------------------------------------------------------------

void MsgInitial::set_dataset
(int n4[4], double h4[4],
 int nx, int ny, int nz,
 int IX, int IY, int IZ)
{
  for (int i=0; i<4; i++) {
    n4_[i] = n4[i];
    h4_[i] = h4[i];
  }
  nx_ = nx;  ny_ = ny;  nz_ = nz;
  IX_ = IX;  IY_ = IY;  IZ_ = IZ;
}

//----------------------------------------------------------------------

void MsgInitial::get_dataset
(int n4[4], double h4[4],
 int * nx, int * ny, int * nz,
 int * IX, int * IY, int * IZ)
{
  for (int i=0; i<4; i++) {
    n4[i] = n4_[i];
    h4[i] = h4_[i];
  }
  (*nx) = nx_;  (*ny) = ny_;  (*nz) = nz_;
  (*IX) = IX_;  (*IY) = IY_;  (*IZ) = IZ_;
}

//----------------------------------------------------------------------

void MsgInitial::set_particle_data
(std::string particle_name, std::string particle_attribute,
 char * data, int data_size, int data_precision)
{
  data_type_      = "particle";
  data_name_      = particle_name;
  data_attribute_ = particle_attribute;

  copy_data_(data,data_size,data_precision);
}

//----------------------------------------------------------------------

void MsgInitial::get_particle_data
(std::string * particle_name,
 std::string * particle_attribute,
 char ** data, int * data_size, int * data_precision)
{
  (*particle_name) = data_name_;
  (*particle_attribute) = data_attribute_;
  (*data) = data_values_;
  const int bytes_per_element = cello::sizeof_precision(data_precision_);
  (*data_size) = data_bytes_ / bytes_per_element;
  (*data_precision) = data_precision_;
}

//======================================================================

void MsgInitial::copy_data_( char * data, int data_size, int data_precision)
{
  // create copy of data
  const int bytes_per_element = cello::sizeof_precision(data_precision);
  data_precision_ = data_precision;
  data_bytes_ = data_size*bytes_per_element;
  data_values_ = new char[data_bytes_];
  data_delete_ = true;
  std::copy_n( data, data_bytes_, data_values_);
}

