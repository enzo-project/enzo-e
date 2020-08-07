// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"

long DataMsg::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

int DataMsg::data_size () const
{
  FieldFace    * ff = field_face_;
  ParticleData * pd = particle_data_;
  auto & fd = face_fluxes_list_;

  bool any_fluxes = fd.size() > 0;

  //--------------------------------------------------
  //  1. determine buffer size (must be consistent with #3)
  //--------------------------------------------------

  Field field (cello::field_descr(), field_data_);

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (ff) ? ff->num_bytes_array(field) : 0;
  const int n_pd = (pd) ? pd->data_size(cello::particle_descr()) : 0;
  const int n_fd = fd.size();

  int size = 0;

  size += sizeof(int); // n_ff_
  size += sizeof(int); // n_fa_
  size += sizeof(int); // n_pd_
  size += sizeof(int); // n_fd_

  size += n_ff;
  size += n_fa;
  size += n_pd;
  for (int i=0; i<n_fd; i++) {
    size += sizeof(int); // face_fluxes_delete_[i] 
    size += fd[i]->data_size();
  }

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

  Field field (cello::field_descr(), field_data_);

  FieldFace    * ff = field_face_;
  ParticleData * pd = particle_data_;
  auto & fd = face_fluxes_list_;

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (ff) ? ff->num_bytes_array(field) : 0;
  const int n_pa = (pd) ? pd->data_size(cello::particle_descr()) : 0;
  const int n_fd = fd.size();

  (*pi++) = n_ff;
  (*pi++) = n_fa;
  (*pi++) = n_pa;
  (*pi++) = n_fd;

  // save field face
  if (n_ff > 0) {
    pc = ff->save_data (pc);
  }
  // save field array
  if (n_ff > 0 && n_fa > 0) {
    ff->face_to_array(field,pc);
    pc += n_fa;
  }
  // save particle data
  if (n_pa > 0) {
    pc = pd->save_data(cello::particle_descr(),pc);
  }
  // save fluxes
  if (n_fd > 0) {
    for (int i=0; i<n_fd; i++) {
      (*pi++) = face_fluxes_delete_[i];
      pc = fd[i]->save_data(pc);
    }
  }

  ASSERT2 ("DataMsg::save_data()",
  	   "Expecting buffer size %d actual size %d",
  	   data_size(),(pc-buffer),
  	   (data_size() == (pc-buffer)));
  
  // return first byte after filled buffer
  return pc;
}

//----------------------------------------------------------------------

char * DataMsg::load_data (char * buffer)
{
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
  const int n_fd = (*pi++);

  // load field face
  if (n_ff > 0) {
    pc = field_face_->load_data (pc);
  }

  // load field array
  if (n_fa > 0) {
    field_array_ = pc;
    pc += n_fa;
  } else {
    field_array_ = NULL;
  }

  // load particle data
  if (n_pa > 0) {
    ParticleData * pd = particle_data_ = new ParticleData;
    pd->allocate(cello::particle_descr());
    pc = pd->load_data(cello::particle_descr(),pc);
  } else {
    particle_data_ = NULL;
  }

  // load flux data
  if (n_fd > 0) {
    face_fluxes_list_.resize(n_fd);
    face_fluxes_delete_.resize(n_fd);
    for (int i=0; i<n_fd; i++) {
      face_fluxes_delete_[i] = (*pi++);
      FaceFluxes * ff = new FaceFluxes;
      face_fluxes_list_[i] = ff;
      pc = ff->load_data(pc);
    }
  }

  return pc;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data, bool is_local)
{
  ParticleData * pd = particle_data_;
  FieldFace    * ff = field_face_;
  char         * fa = field_array_;

  // Update particles
  
  if (pd != NULL) {

    // Insert new particles 

    Particle particle = data->particle();

    int count = 0;
    for (int it=0; it<particle.num_types(); it++) {
      count += particle.gather (it, 1, &pd);
    }
    cello::simulation()->data_insert_particles(count);
    
    delete particle_data_;
    particle_data_ = NULL; 
  }
  
  // Update fields
  
  if (ff != NULL && fa != NULL) {

    Field field_dst = data->field();

    if (is_local) {

      Field field_src(cello::field_descr(),field_data_);

      ff->face_to_face(field_src, field_dst);

    } else { // ! is_local

      // invert face since incoming not outgoing

      ff->invert_face();

      ff->array_to_face(fa,field_dst);

    }
  }

  // Update fluxes

  FluxData * flux_data = data->flux_data();
  if (face_fluxes_list_.size() > 0) {
    for (int i=0; i<face_fluxes_list_.size(); i++) {
      FaceFluxes * face_fluxes = face_fluxes_list_[i];
      Face face = face_fluxes->face();
      flux_data->sum_neighbor_fluxes
        (face_fluxes,face.axis(), 1 - face.face(), i);
    }
    face_fluxes_list_.clear();
  }
  if (ff != NULL) {
    delete field_face_;
    field_face_ = NULL;
  }

}
