// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"

long DataMsg::counter[CONFIG_NODE_SIZE] = {0};

#define CHECK

// #define TRACE_COARSE_ARRAY_VALUES

//----------------------------------------------------------------------

void DataMsg::set_coarse_array
  (Field field,
   int iam3[3], int iap3[3],
   int ifms3[3], int ifps3[3],
   int ifmr3[3], int ifpr3[3],
   const std::vector<int> & coarse_field_list)
{
  for (int i=0; i<3; i++) {
    iam3_cf_[i] = iam3[i];
    iap3_cf_[i] = iap3[i];
    ifms3_cf_[i] = ifms3[i];
    ifps3_cf_[i] = ifps3[i];
    ifmr3_cf_[i] = ifmr3[i];
    ifpr3_cf_[i] = ifpr3[i];
  }
  coarse_field_list_ = coarse_field_list;

  const int nf3[3] = {(ifps3[0] - ifms3[0]),
                      (ifps3[1] - ifms3[1]),
                      (ifps3[2] - ifms3[2])};
  const int na3[3] = {(iap3[0] - iam3[0]),
                      (iap3[1] - iam3[1]),
                      (iap3[2] - iam3[2])};
  const int nf = coarse_field_list.size();
  const int na = na3[0]*na3[1]*na3[2];

  coarse_field_buffer_.resize(nf*na);
  std::fill(coarse_field_buffer_.begin(),coarse_field_buffer_.end(),0.0);
  for (int i_f=0; i_f<nf; i_f++) {
    const int index_field = coarse_field_list[i_f];
    cello_float * coarse_field = coarse_field_buffer_.data() + i_f*na;
    // get Field dimensions
    int mfx,mfy,mfz;
    field.dimensions(index_field,&mfx,&mfy,&mfz);
    // get field values
    cello_float * field_values = (cello_float *)field.values(index_field);
    // determine starting offset for field
    const int if0 = ifms3[0] + mfx*(ifms3[1] + mfy*ifms3[2]);
    const cello_float rr = 1.0;
    // compute cell width ratio r
    const float r = nf3[0] / na3[0];
#ifdef CHECK
    ASSERT1 ("DataMsg::set_coarse_array",
             "Field-to-coarse array axis ratio r=%g is not 1.0 or 2.0",
             r, (r==1.0 || r==2.0));
#endif    
    // compute constant factor v for r!=1
    const float v = (r==1)?1.0 : 1.0/cello::num_children();
    for (int kz=0; kz<nf3[2]; kz++) {
      for (int ky=0; ky<nf3[1]; ky++) {
        for (int kx=0; kx<nf3[0]; kx++) {
          int ka = (kx/r) + na3[0]*((ky/r) + na3[1]*(kz/r));
          int kf = if0 + kx + mfx*(ky + mfy*kz);
          int kf2 = (ifms3[0] + kx) + mfx*((ifms3[1]+ky) + mfy*(ifms3[2]+kz));
          coarse_field[ka] += rr*field_values[kf];
        }
      }
    }
#ifdef TRACE_COARSE_ARRAY_VALUES
    if (i_f == 0) {
      CkPrintf ("COARSE_ARRAY_VALUES set_coarse_array  field %d (%d %d %d)\n",
                i_f,na3[0],na3[1],na3[2]);
      for (int iz=0; iz<na3[2]; iz++) {
        for (int iy=0; iy<na3[1]; iy++) {
          CkPrintf ("COARSE_ARRAY_VALUES coarse_field %d %d %d: ",0,iy,iz);
          for (int ix=0; ix<na3[0]; ix++) {
            int i = ix+ na3[0]*(iy+ na3[1]*iz);
            CkPrintf (" %6.3g",coarse_field[i]);
          }
          CkPrintf ("\n");
        }
      }
    }
#endif
  }
}

//----------------------------------------------------------------------

int DataMsg::data_size () const
{
  FieldFace    * ff = field_face_;
  ParticleData * pd = particle_data_;
  auto & fd = face_fluxes_list_;

  //--------------------------------------------------
  //  1. determine buffer size (must be consistent with #3)
  //--------------------------------------------------

  int debug_counter = 0;
  
  Field field (cello::field_descr(), field_data_u_);

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (ff) ? ff->num_bytes_array(field) : 0;
  const int n_pd = (pd) ? pd->data_size(cello::particle_descr()) : 0;
  const int n_fd = fd.size();

  int size = 0;

  SIZE_SCALAR_TYPE(size,int,n_ff);
  SIZE_SCALAR_TYPE(size,int,n_fa);
  SIZE_SCALAR_TYPE(size,int,n_pd);
  SIZE_SCALAR_TYPE(size,int,n_fd);

  size += n_ff;
  size += n_fa;
  size += n_pd;

  if (n_fd > 0) {
    size += n_fd*sizeof(int); // face_fluxes_delete_[i]
    for (int i=0; i<n_fd; i++) {
      size += fd[i]->data_size();
    }
  }

  // Coarse array for interpolation

  const int nax=iap3_cf_[0]-iam3_cf_[0];
  const int nay=iap3_cf_[1]-iam3_cf_[1];
  const int naz=iap3_cf_[2]-iam3_cf_[2];
  const int na = nax*nay*naz;

  SIZE_SCALAR_TYPE(size,int,na);  // coarse array size (0 if none)

  if (na > 0) {

    SIZE_VECTOR_TYPE(size,int,coarse_field_list_);
    SIZE_VECTOR_TYPE(size,cello_float,coarse_field_buffer_);

    size += 3*sizeof(int); // iam3_cf_
    size += 3*sizeof(int); // iap3_cf_
    size += 3*sizeof(int); // ifms3_cf_
    size += 3*sizeof(int); // ifps3_cf_
    size += 3*sizeof(int); // ifmr3_cf_
    size += 3*sizeof(int); // ifpr3_cf_
  }
  size += 8*sizeof(char); // tag_
  
  return size;
}

//----------------------------------------------------------------------

char * DataMsg::save_data (char * buffer) const
{
  union {
    char * pc;
    int  * pi;
    cello_float * pcf;
  };

  pc = buffer;
  int debug_counter = 0;

  Field field (cello::field_descr(), field_data_u_);

  FieldFace    * ff = field_face_;
  ParticleData * pd = particle_data_;
  auto & fd = face_fluxes_list_;

  const int n_ff = (ff) ? ff->data_size() : 0;
  const int n_fa = (ff) ? ff->num_bytes_array(field) : 0;
  const int n_pa = (pd) ? pd->data_size(cello::particle_descr()) : 0;
  const int n_fd = fd.size();

  SAVE_SCALAR_TYPE(pc,int,n_ff);
  SAVE_SCALAR_TYPE(pc,int,n_fa);
  SAVE_SCALAR_TYPE(pc,int,n_pa);
  SAVE_SCALAR_TYPE(pc,int,n_fd);

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

  // Coarse face array for interpolation

  const int nax=iap3_cf_[0]-iam3_cf_[0];
  const int nay=iap3_cf_[1]-iam3_cf_[1];
  const int naz=iap3_cf_[2]-iam3_cf_[2];
  const int na = nax*nay*naz;

  SAVE_SCALAR_TYPE(pc,int,na);

  if (na > 0) {

    SAVE_VECTOR_TYPE(pc,int,coarse_field_list_);
    SAVE_VECTOR_TYPE(pc,cello_float,coarse_field_buffer_);

    for (int i=0; i<3; i++) {
      (*pi++) = iam3_cf_[i];
      (*pi++) = iap3_cf_[i];
      (*pi++) = ifms3_cf_[i];
      (*pi++) = ifps3_cf_[i];
      (*pi++) = ifmr3_cf_[i];
      (*pi++) = ifpr3_cf_[i];
    }
  
  }
  strncpy (pc,tag_,8);
  pc+=8*sizeof(char);

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
    cello_float * pcf;
  };

  pc = buffer;
  int debug_counter = 0;

  int n_ff,n_fa,n_pa,n_fd;
  LOAD_SCALAR_TYPE(pc,int,n_ff);
  LOAD_SCALAR_TYPE(pc,int,n_fa);
  LOAD_SCALAR_TYPE(pc,int,n_pa);
  LOAD_SCALAR_TYPE(pc,int,n_fd);

  // load field face
  if (n_ff > 0) {
    field_face_ = new FieldFace(cello::rank());
    pc = field_face_->load_data (pc);
  } else {
    field_face_ = nullptr;
  }

  // load field array
  if (n_fa > 0) {
    field_array_u_ = pc;
    pc += n_fa;
  } else {
    field_array_u_ = nullptr;
  }

  // load particle data
  if (n_pa > 0) {
    ParticleData * pd = particle_data_ = new ParticleData;
    pd->allocate(cello::particle_descr());
    pc = pd->load_data(cello::particle_descr(),pc);
  } else {
    particle_data_ = nullptr;
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

  int na;
  LOAD_SCALAR_TYPE(pc,int,na);

  if (na > 0) {

    LOAD_VECTOR_TYPE(pc,int,coarse_field_list_);
    LOAD_VECTOR_TYPE(pc,cello_float,coarse_field_buffer_);

    for (int i=0; i<3; i++) {
      iam3_cf_[i] = (*pi++);
      iap3_cf_[i] = (*pi++);
      ifms3_cf_[i] = (*pi++);
      ifps3_cf_[i] = (*pi++);
      ifmr3_cf_[i] = (*pi++);
      ifpr3_cf_[i] = (*pi++);
    }
  }
  strncpy (tag_,pc,8);
  tag_[8] = 0;
  pc+=8*sizeof(char);
    
  return pc;
}

//----------------------------------------------------------------------

void DataMsg::update (Data * data, bool is_local)
{
  ParticleData * pd = particle_data_;
  FieldFace    * ff = field_face_;
  char         * fa = field_array_u_;

  // Update particles
  if (pd != nullptr) {

    // Insert new particles 

    Particle particle = data->particle();

    int count = 0;
    for (int it=0; it<particle.num_types(); it++) {
      count += particle.gather (it, 1, &pd);
    }
    cello::simulation()->data_insert_particles(count);
    
    delete particle_data_;
    particle_data_ = nullptr; 
  }
  
  // Update fields
  
  if (ff != nullptr && fa != nullptr) {

    Field field_dst = data->field();

    if (is_local) {

      Field field_src(cello::field_descr(),field_data_u_);

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
  if (ff != nullptr) {
    delete field_face_;
    field_face_ = nullptr;
  }

  // Updated coarse array
  const int na3[3] =
    {(iap3_cf_[0] - iam3_cf_[0]),
     (iap3_cf_[1] - iam3_cf_[1]),
     (iap3_cf_[2] - iam3_cf_[2])};
  
  const int na = na3[0]*na3[1]*na3[2];

  if (na>0) {

    const int nf3[3] =
      {(ifps3_cf_[0] - ifms3_cf_[0]),
       (ifps3_cf_[1] - ifms3_cf_[1]),
       (ifps3_cf_[2] - ifms3_cf_[2])};
    const int na3[3] =
      {(iap3_cf_[0] - iam3_cf_[0]),
       (iap3_cf_[1] - iam3_cf_[1]),
       (iap3_cf_[2] - iam3_cf_[2])};
    
    // const int ia_start = iaxm + mx*(iaym + my*iazm);
    const int nf = coarse_field_list_.size();
    
    for (int i_f=0; i_f<nf; i_f++) {

      Field field = data->field();
      const int index_field = coarse_field_list_[i_f];
      int m3_c[3];
      field.coarse_dimensions(index_field,m3_c,m3_c+1,m3_c+2);
      const int ic0 = iam3_cf_[0] + m3_c[0]*(iam3_cf_[1] + m3_c[1]*iam3_cf_[2]);

      cello_float * coarse_buffer = coarse_field_buffer_.data() + i_f*na;
      cello_float * coarse_field = (cello_float *)field.coarse_values(index_field);
      
      for (int kz=0; kz<nf3[2]; kz++) {
        for (int ky=0; ky<nf3[1]; ky++) {
          for (int kx=0; kx<nf3[0]; kx++) {
            const int kb = kx + nf3[0]*(ky +  nf3[1]*kz);
            const int kc = ic0 + kx + m3_c[0]*(ky + m3_c[1]*kz);
            coarse_field[kc] = coarse_buffer[kb];
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void DataMsg::print (const char * message) const
{
  CkPrintf ("%s %d --------------------------------------------------\n",
            message,tag_);
  CkPrintf ("%s %s DATA_MSG %p\n",  message,tag_,(void *)this);
  CkPrintf ("%s %s DATA_MSG field_face_    = %p\n",
            message,tag_,(void*)field_face_);
  CkPrintf ("%s %s DATA_MSG field_data_u_    = %p\n",
            message,tag_,(void*)field_data_u_);
  CkPrintf ("%s %s DATA_MSG particle_data_ = %p\n",
            message,tag_,(void*)particle_data_);
  CkPrintf ("%s %s DATA_MSG particle_data_delete_ = %d\n",
            message,tag_,particle_data_delete_?1:0);
  CkPrintf ("%s %s DATA_MSG |face_fluxes_list_| = %lu\n",
            message,tag_,face_fluxes_list_.size());
  CkPrintf ("%s %s DATA_MSG |face_fluxes_delete_| = %lu\n",
            message,tag_,face_fluxes_delete_.size());
  CkPrintf ("%s %s DATA_MSG field_face_delete_ = %d\n",
            message,tag_,field_face_delete_?1:0);
  CkPrintf ("%s %s DATA_MSG field_data_delete_ = %d\n",
            message,tag_,field_data_delete_?1:0);
  CkPrintf ("%s %s DATA_MSG coarse_field_buffer_.sum = %f\n", message,tag_,
            std::accumulate(coarse_field_buffer_.begin(),coarse_field_buffer_.end(),0.0));
  CkPrintf ("%s %s DATA_MSG coarse_field_list_.sum = %d\n", message,tag_,
            std::accumulate
            (coarse_field_list_.begin(),
             coarse_field_list_.end(),0));
  CkPrintf ("%s %s DATA_MSG coarse iam3_cf_   = %d %d %d\n",
            message,tag_,iam3_cf_[0],iam3_cf_[1],iam3_cf_[2]);
  CkPrintf ("%s %s DATA_MSG coarse iap3_cf_   = %d %d %d\n",
            message,tag_,iap3_cf_[0],iap3_cf_[1],iap3_cf_[2]);
  CkPrintf ("%s %s DATA_MSG coarse ifms3_cf_   = %d %d %d\n",
            message,tag_,ifms3_cf_[0],ifms3_cf_[1],ifms3_cf_[2]);
  CkPrintf ("%s %s DATA_MSG coarse ifps3_cf_   = %d %d %d\n",
            message,tag_,ifps3_cf_[0],ifps3_cf_[1],ifps3_cf_[2]);
  CkPrintf ("%s %s DATA_MSG coarse ifmr3_cf_   = %d %d %d\n",
            message,tag_,ifmr3_cf_[0],ifmr3_cf_[1],ifmr3_cf_[2]);
  CkPrintf ("%s %s DATA_MSG coarse ifpr3_cf_   = %d %d %d\n",
            message,tag_,ifpr3_cf_[0],ifpr3_cf_[1],ifpr3_cf_[2]);
    
  fflush(stdout);
}
