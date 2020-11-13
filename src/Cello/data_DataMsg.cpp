// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"

long DataMsg::counter[CONFIG_NODE_SIZE] = {0};

// #define DEBUG_PADDED_ARRAY
//#define ACTIVATE_PADDED_ARRAY
#define CHECK

//----------------------------------------------------------------------

void DataMsg::set_padded_face
(int if3[3],int ma3[3],
 int iam3[3], int iap3[3],
 int ifm3[3], int ifp3[3],
 std::vector<int> padded_face_field_list, Field field)
{
#ifdef ACTIVATE_PADDED_ARRAY  
  for (int i=0; i<3; i++) {
    if3_pf_[i] = if3[i];
    m3_pf_[i] = ma3[i];
    im3_pf_[i] = iam3[i];
    ip3_pf_[i] = iap3[i];
  }
  ASSERT6("DataMsg::set_padded_face",
          "field section vs array size discrepancy: field (%d %d %d) != array (%d %d %d)",
          (ifp3[0]-ifm3[0]),(ifp3[1]-ifm3[1]),(ifp3[2]-ifm3[2]),
          (iap3[0]-iam3[0]),(iap3[1]-iam3[1]),(iap3[2]-iam3[2]),
          ((ifp3[0]-ifm3[0]) == (iap3[0] - iam3[0])) &&
          ((ifp3[1]-ifm3[1]) == (iap3[1] - iam3[1])) &&
          ((ifp3[2]-ifm3[2]) == (iap3[2] - iam3[2])));
          
  const int nx=(ifp3[0]-ifm3[0]);
  const int ny=(ifp3[1]-ifm3[1]);
  const int nz=(ifp3[2]-ifm3[2]);
  
  padded_face_field_list_ = padded_face_field_list;
  const int ma = ma3[0]*ma3[1]*ma3[2];
  const int nf = padded_face_field_list.size();
  padded_face_.resize(nf*ma);
  for (int i_f=0; i_f<nf; i_f++) {
    const int index_field = padded_face_field_list[i_f];
    int mf3[3];
    field.dimensions(index_field,mf3,mf3+1,mf3+2);
    const int mf = mf3[0]*mf3[1]*mf3[2];
    auto field_values = (cello_float *)field.values(index_field);
    int ia0 = iam3[0] + ma3[0]*(iam3[1] + ma3[1]*iam3[2]);
    int if0 = ifm3[0] + mf3[0]*(ifm3[1] + mf3[1]*ifm3[2]);
    for (int kz=0; kz<nz; kz++) {
      for (int ky=0; ky<ny; ky++) {
        for (int kx=0; kx<nx; kx++) {
          int ka=kx + ma3[0]*(ky + ma3[1]*kz) + ia0;
          int kf=kx + mf3[0]*(ky + mf3[1]*kz) + if0;
#ifdef CHECK
          ASSERT6 ("DataMsg::set_padded_face",
                  "padded_face_ array index (%d %d %d) out of range (%d %d %d)",
                   kx+iam3[0],ky+iam3[1],kz+iam3[2],
                   ma3[0],ma3[1],ma3[2],
                   (0 <= ka && ka < ma));
          ASSERT6 ("DataMsg::set_padded_face",
                  "Field index (%d %d %d) out of range (%d %d %d)",
                   kx+ifm3[0],ky+ifm3[1],kz+ifm3[2],
                   mf3[0],mf3[1],mf3[2],
                   (0 <= ka && ka < mf));
#endif          
          padded_face_[ka] = field_values[kf];
        }
      }
    }
  }
#endif  
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

  // Padded array for interpolation
  size += sizeof(int); // padded array size (0 if none)
  const int m = (m3_pf_[0]*m3_pf_[1]*m3_pf_[2]);
  if (m > 0) {
    size += sizeof(int); // padded_face_field_list_.size()

    int n = padded_face_field_list_.size();
    size += n*sizeof(int); // padded_face_field_list_
    size += (n*m)*sizeof(cello_float); // padded_face_

    size += 3*sizeof(int); // m3_pf_[3]
    size += 3*sizeof(int); // if3_pf_[3]
    size += 3*sizeof(int); // im3_pf_[3]
    size += 3*sizeof(int); // ip3_pf_[3]
  }
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("DEBUG DataMsg::data_size m3 %d %d %d\n",
            m3_pf_[0],m3_pf_[1],m3_pf_[2]);
#endif
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

  // Padded face array for interpolation

  const int m = (m3_pf_[0]*m3_pf_[1]*m3_pf_[2]);
  (*pi++) = m;
  if (m > 0) {
    const int n = padded_face_field_list_.size();
    (*pi++) = n;
    for (int i=0; i<n; i++) {
      (*pi++) = padded_face_field_list_[i];
    }
    for (int i=0; i<(n*m); i++) {
      (*pcf++) = padded_face_[i];
    }

    for (int i=0; i<3; i++) {
      (*pi++) = m3_pf_[i];
      (*pi++) = if3_pf_[i];
      (*pi++) = im3_pf_[i];
      (*pi++) = ip3_pf_[i];
    }
#ifdef DEBUG_PADDED_ARRAY
    CkPrintf ("DEBUG DataMsg::save_data m3 %d %d %d\n",
              m3_pf_[0],m3_pf_[1],m3_pf_[2]);
#endif
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
    cello_float * pcf;
  };

  pc = buffer;

  field_face_ = new FieldFace(cello::rank());

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

  const int m = (*pi++);
  if (m > 0) {
    const int n = (*pi++);
    padded_face_field_list_.resize(n);
    for (int i=0; i<n; i++) {
      padded_face_field_list_[i] = (*pi++);
    }
    padded_face_.resize(n*m);
    for (int i=0; i<(n*m); i++) {
      padded_face_[i] = (*pcf++);
    }

    for (int i=0; i<3; i++) {
      m3_pf_[i] = (*pi++);
      if3_pf_[i] = (*pi++);
      im3_pf_[i] = (*pi++);
      ip3_pf_[i] = (*pi++);
    }
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG DataMsg::load_data m3 %d %d %d\n",
              m3_pf_[0],m3_pf_[1],m3_pf_[2]);
#endif
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

  // Updated padded array

  const int m = m3_pf_[0]*m3_pf_[1]*m3_pf_[2];
  if (m>0) {
    const int mx(m3_pf_[0]), my(m3_pf_[1]), mz(m3_pf_[2]);
    const int ixm(im3_pf_[0]), iym(im3_pf_[1]), izm(im3_pf_[2]);
    const int ixp(ip3_pf_[0]), iyp(ip3_pf_[1]), izp(ip3_pf_[2]);
    const int ifx(if3_pf_[0]), ify(if3_pf_[1]), ifz(if3_pf_[2]);

    auto array = data->padded_face_array_allocate (ifx,ify,ifz, mx,my,mz);
    
    for (int iz=izm; iz<izp; iz++) {
      for (int iy=iym; iy<iyp; iy++) {
        for (int ix=ixm; ix<ixp; ix++) {
          int i = ix + mx*(iy + my*iz);
          array[i] = padded_face_[i];
        }
      }
    }
#ifdef DEBUG_PADDED_ARRAY
    CkPrintf ("DEBUG DataMsg::update m3 %d %d %d\n",
              m3_pf_[0],m3_pf_[1],m3_pf_[2]);
#endif
  }

}
