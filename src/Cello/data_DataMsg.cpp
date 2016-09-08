// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_NEW_REFRESH

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

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  ParticleDescr * particle_descr = simulation->particle_descr();

  Field field (field_descr, field_data_);

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (fa) ? ff->num_bytes_array(field) : 0;
  const int n_pa = (pd) ? pd->data_size(particle_descr) : 0;

  int size = 0;

  size += sizeof(int); // n_ff_
  size += sizeof(int); // n_fa_
  size += sizeof(int); // n_pa_

  size += n_ff*sizeof(char);
  size += n_fa*sizeof(char);
  size += n_pa*sizeof(char);

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg data_size %d\n",this,size);
#endif
  return size;
}

//----------------------------------------------------------------------

char * DataMsg::save_data (char * buffer) const
{

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  ParticleDescr * particle_descr = simulation->particle_descr();

  Field field (field_descr, field_data_);

  FieldFace    * ff = field_face_;
  char         * fa = field_array_;
  ParticleData * pd = particle_data_;

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (fa) ? ff->num_bytes_array(field) : 0;
  const int n_pa = (pd) ? pd->data_size(particle_descr) : 0;

  (*pi++) = n_ff;
  (*pi++) = n_fa;
  (*pi++) = n_pa;

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg save_data n_ff %d\n",this, n_ff);
#endif
  if (n_ff > 0) {
    pc = ff->save_data (pc);
  }
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg save_data n_fa %d\n",this, n_fa);
#endif
  if (n_ff > 0 && n_fa > 0) {
    ff->face_to_array(field,pc);
    pc += n_fa;
  }
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg save_data n_pa %d\n",this, n_pa);
#endif
  if (n_pa > 0) {
    pc = pd->save_data(particle_descr,pc);
  }

  // return first byte after filled buffer

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg save_data %d\n",this, (pc-buffer));
#endif
  return pc;
}

//----------------------------------------------------------------------

char * DataMsg::load_data (char * buffer)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  ParticleDescr * particle_descr = simulation->particle_descr();

  // 2. De-serialize message data from input buffer into the allocated
  // message (must be consistent with pack())
 
  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  field_face_ = new FieldFace;

  const int n_ff = (*pi++);
  const int n_fa = (*pi++);
  const int n_pa = (*pi++);

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg load_data n_ff %d\n",this, n_ff);
  fflush(stdout);
#endif
  if (n_ff > 0) {
    pc = field_face_->load_data (pc);
  }

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg load_data n_fa %d\n",this, n_fa);
  fflush(stdout);
#endif
  if (n_fa > 0) {
    field_array_ = pc;
    pc += n_fa;
  } else {
    field_array_ = NULL;
  }

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg load_data n_pa %d\n",this, n_pa);
  fflush(stdout);
#endif
  if (n_pa > 0) {
    ParticleData * pd = particle_data_ = new ParticleData;

    pd->allocate(particle_descr);
    pc = pd->load_data(particle_descr,pc);
  } else {
    particle_data_ = NULL;
  }

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%p DataMsg load_data %d\n",this, (pc-buffer));
#endif
  return pc;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data, bool is_local)
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr    *    field_descr = simulation->   field_descr();
 
  Field field_dst = data->field();
 
  FieldData    * fd = field_data();
  ParticleData * pd = particle_data();
  FieldFace    * ff = field_face();
  char         * fa = field_array();

  if (pd != NULL) {

    // Insert new particles 

    Particle particle = data->particle();

    int count = 0;
    for (int it=0; it<particle.num_types(); it++) {
      count += particle.gather (it, 1, &pd);
    }
    simulation->monitor_insert_particles(count);
    
    delete_particle_data();
  }

  if (fa != NULL) {

    if (is_local) {

      Field field_src(field_descr,fd);
      ff->face_to_face(field_src, field_dst);

    } else { // ! is_local

      // Invert face since incoming not outgoing

      ff->invert_face();

      ff->array_to_face(fa,field_dst);

    }
    delete field_face_;
    field_face_ = NULL;

  }
}
