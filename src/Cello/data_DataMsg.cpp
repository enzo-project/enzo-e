// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"

// #define DEBUG_DATA_MSG

long DataMsg::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

int DataMsg::data_size () const
{
  FieldFace    * ff = field_face_;
  char         * fa = field_array_;
  ParticleData * pd = particle_data_;

  //--------------------------------------------------
  //  1. determine buffer size (must be consistent with #3)
  //--------------------------------------------------

  FieldDescr    * field_descr    = cello::field_descr();
  ParticleDescr * particle_descr = cello::particle_descr();

  Field field (field_descr, field_data_);

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (ff) ? ff->num_bytes_array(field) : 0;
  const int n_pa = (pd) ? pd->data_size(particle_descr) : 0;

  int size = 0;

  size += sizeof(int); // n_ff_
  size += sizeof(int); // n_fa_
  size += sizeof(int); // n_pa_

  size += n_ff*sizeof(char);
  size += n_fa*sizeof(char);
  size += n_pa*sizeof(char);

  return size;
}

//----------------------------------------------------------------------

char * DataMsg::save_data (char * buffer) const
{

#ifdef DEBUG_DATA_MSG
  CkPrintf ("%d %s:%d DEBUG_DATA_MSG save_data %p\n",
	    CkMyPe(),__FILE__,__LINE__,this);
  fflush(stdout);
#endif

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  FieldDescr    * field_descr    = cello::field_descr();
  ParticleDescr * particle_descr = cello::particle_descr();

  Field field (field_descr, field_data_);

  FieldFace    * ff = field_face_;
  char         * fa = field_array_;
  ParticleData * pd = particle_data_;

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (ff) ? ff->num_bytes_array(field) : 0;
  const int n_pa = (pd) ? pd->data_size(particle_descr) : 0;

  (*pi++) = n_ff;
  (*pi++) = n_fa;
  (*pi++) = n_pa;

  if (n_ff > 0) {
    pc = ff->save_data (pc);
  }
  if (n_ff > 0 && n_fa > 0) {
    ff->face_to_array(field,pc);
    pc += n_fa;
  }
  if (n_pa > 0) {
    pc = pd->save_data(particle_descr,pc);
  }

  // return first byte after filled buffer


  //  ASSERT2 ("DataMsg::save_data()",
  //	   "Expecting buffer size %d actual size %d",
  //	   data_size(),(pc-buffer),
  //	   (data_size() == (pc-buffer)));
  
  return pc;
}

//----------------------------------------------------------------------

char * DataMsg::load_data (char * buffer)
{
#ifdef DEBUG_DATA_MSG
  CkPrintf ("%d %s:%d DEBUG_DATA_MSG load_data %p\n",
	    CkMyPe(),__FILE__,__LINE__,this);
  fflush(stdout);
#endif

  ParticleDescr * particle_descr = cello::particle_descr();

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())
 
  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  field_face_ = new FieldFace;
#ifdef DEBUG_FIELD_FACE
  CkPrintf ("%d %s:%d DEBUG_FIELD_FACE creating %p\n",CkMyPe(),__FILE__,__LINE__,field_face_);
#endif  

  const int n_ff = (*pi++);
  const int n_fa = (*pi++);
  const int n_pa = (*pi++);

  if (n_ff > 0) {
    pc = field_face_->load_data (pc);
  }

  if (n_fa > 0) {
    field_array_ = pc;
    pc += n_fa;
  } else {
    field_array_ = NULL;
  }

  if (n_pa > 0) {
    ParticleData * pd = particle_data_ = new ParticleData;
    pd->allocate(particle_descr);
    pc = pd->load_data(particle_descr,pc);
  } else {
    particle_data_ = NULL;
  }

  //  ASSERT2 ("DataMsg::load_data()",
  //  	   "Expecting buffer size %d actual size %d",
  //  	   data_size(),(pc-buffer),
  //  	   (data_size() == (pc-buffer)));

  return pc;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data, bool is_local)
{
  Simulation * simulation  = cello::simulation();
  FieldDescr * field_descr = cello::field_descr();
 
  Field field_dst = data->field();
 
  FieldData    * fd = field_data();
  ParticleData * pd = particle_data();
  FieldFace    * ff = field_face();
  char         * fa = field_array();

#ifdef DEBUG_DATA_MSG
  CkPrintf ("%d %s:%d DEBUG_DATA_MSG update %p (fd:%p pd:%p ff:%p fa:%p local:%d)\n",
	    CkMyPe(),__FILE__,__LINE__,this,fd,pd,ff,fa,is_local);
  fflush(stdout);
#endif


  if (pd != NULL) {

    // Insert new particles 

    Particle particle = data->particle();

    int count = 0;
    for (int it=0; it<particle.num_types(); it++) {
      count += particle.gather (it, 1, &pd);
    }
    simulation->data_insert_particles(count);
    
    delete particle_data_;
    particle_data_ = NULL; 


  }
  

  if (ff != NULL && fa != NULL) {

    if (is_local) {

      Field field_src(field_descr,fd);

      ff->face_to_face(field_src, field_dst);

    } else { // ! is_local

      // Invert face since incoming not outgoing

      ff->invert_face();

      ff->array_to_face(fa,field_dst);

    }
  }
#ifdef DEBUG_DATA_MSG
  CkPrintf ("%d %s:%d DEBUG_DATA_MSG delete? %p (%p %d)\n",
	    CkMyPe(),__FILE__,__LINE__,this,ff,field_face_delete_);
  fflush(stdout);
#endif
  if (ff != NULL) {
    delete field_face_;
    field_face_ = NULL;
  }

}
