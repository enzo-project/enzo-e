// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-22
/// @brief    

#include "data.hpp"

long DataMsg::counter[CONFIG_NODE_SIZE] = {0};

// #define DEBUG_PADDED_ARRAY
// #define DEBUG_TRACE_SUMS

#define ACTIVATE_PADDED_ARRAY
// #define TRACE_PADDED_ARRAY_VALUES
#define CHECK

//----------------------------------------------------------------------

void DataMsg::set_padded_face
(int ma3[3],int n3[3],int r, int v,
 int iam3[3], int ifm3[3], int if3[3],
 std::vector<int> padded_face_field_list, Field field,
 std::string block_name)
{
#ifdef ACTIVATE_PADDED_ARRAY  
  for (int i=0; i<3; i++) {
    ma3_pf_[i]  = ma3[i];
    n3_pf_[i]   = n3[i];
    iam3_pf_[i] = iam3[i];
    if3_pf_[i]  = if3[i];
  }
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG DataMsg::set_padded_face %p ma3 %d %d %d\n",
              (void *)this,ma3_pf_[0],ma3_pf_[1],ma3_pf_[2]);
    CkPrintf ("DEBUG DataMsg::set_padded_face %p n3 %d %d %d\n",
              (void *)this,n3_pf_[0],n3_pf_[1],n3_pf_[2]);
    CkPrintf ("DEBUG DataMsg::set_padded_face %p iam3 %d %d %d\n",
              (void *)this,iam3_pf_[0],iam3_pf_[1],iam3_pf_[2]);
    CkPrintf ("DEBUG DataMsg::set_padded_face %p if3 %d %d %d\n",
              (void *)this,if3_pf_[0],if3_pf_[1],if3_pf_[2]);
#endif
  padded_face_field_list_ = padded_face_field_list;

  const int nx(n3_pf_[0]), ny(n3_pf_[1]), nz(n3_pf_[2]);
  const int nf = padded_face_field_list.size();
  const int na(nx*ny*nz);
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE nf %s:%d %d\n",__FILE__,__LINE__,nf);
  CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE na %s:%d %d\n",__FILE__,__LINE__,na);
#endif  

  padded_face_.resize(nf*na);
  std::fill(padded_face_.begin(),padded_face_.end(),0.0);
  cello_float * padded_field = padded_face_.data();
  for (int i_f=0; i_f<nf; i_f++) {
    const int index_field = padded_face_field_list[i_f];
    const int ia0 = i_f*na;
    // get Field dimensions
    int mxf,myf,mzf;
    field.dimensions(index_field,&mxf,&myf,&mzf);
    const int mf = mxf*myf*mzf;
    // get field values
    auto field_values = (cello_float *)field.values(index_field);
    // determine starting offset for array and field
    int if0 = ifm3[0] + mxf*(ifm3[1] + myf*ifm3[2]);
    // copy field values to array, looping over field elements
    const int rank = cello::rank();
    const int nxp = (rank >= 1) ? r*nx : 1;
    const int nyp = (rank >= 2) ? r*ny : 1;
    const int nzp = (rank >= 3) ? r*nz : 1;
    //    const cello_float rr = 1.0/cello::num_children();
    const cello_float rr = 1.0/v;
    for (int kz=0; kz<nzp; kz++) {
      for (int ky=0; ky<nyp; ky++) {
        for (int kx=0; kx<nxp; kx++) {
          int ka=ia0 + (kx/r) + nx*((ky/r) + ny*(kz/r));
          int kf=if0 + kx + mxf*(ky + myf*kz);
#ifdef CHECK
          ASSERT8 ("DataMsg::set_padded_face",
                  "padded_face_ array index (%d %d %d) out of range (%d %d %d) %d %d",
                   kx/r,ky/r,kz/r,
                   nx,ny,nz,ka-ia0,na,
                   (0 <= ka && ka < na*nf));
          ASSERT6 ("DataMsg::set_padded_face",
                  "Field index (%d %d %d) out of range (%d %d %d)",
                   kx+ifm3[0],ky+ifm3[1],kz+ifm3[2],
                   mxf,myf,mzf,
                   (0 <= kf && kf < mf));
#endif
          padded_field[ka] += rr*field_values[kf];
        }
      }
    }
#ifdef DEBUG_TRACE_SUMS
    {
      int count = 0;
      double sum = 0.0;
      int zero = 0;
      for (int iz=0; iz<nz; iz++) {
        for (int iy=0; iy<ny; iy++) {
          for (int ix=0; ix<nx; ix++) {
            int i=ia0+ix+nx*(iy+ny*iz);
            count++;
            sum += padded_face_[i];
            if (padded_face_[i] == 0) ++zero;
          }
        }
      }
      CkPrintf ("DEBUG_TRACE_SUM %s send buffer [%d %d %d] n3 %d %d %d i_f %d %s %d %d %8.4g\n",
                block_name.c_str(), if3[0],if3[1],if3[2],
                nx,ny,nz,i_f,tag_,zero,count,sum);
    }
#endif    
#ifdef TRACE_PADDED_ARRAY_VALUES
    if (i_f == 0) {
      CkPrintf ("PADDED_ARRAY_VALUES %s:%d i_f %d set_padded_array field %p r %d\n",
                __FILE__,__LINE__,i_f,(void*)&field_values[if0],r);
      for (int iz=0; iz<nzp; iz++) {
        for (int iy=0; iy<nyp; iy++) {
          CkPrintf ("PADDED_ARRAY_VALUES set_padded_array field i_f %d %p %d %d %d: ",
                    i_f,(void*)&field_values[if0],0,iy,iz);
          for (int ix=0; ix<nxp; ix++) {
            int i=if0 + ix + mxf*(iy + myf*iz);
            CkPrintf (" %6.3g",field_values[i]);
          }
          CkPrintf ("\n");
        }
      }
      CkPrintf ("PADDED_ARRAY_VALUES %s:%d i_f %d set_padded_array buffer %p\n",
                __FILE__,__LINE__,i_f,(void*)&padded_face_[ia0]);
      for (int iz=0; iz<nz; iz++) {
        for (int iy=0; iy<ny; iy++) {
          CkPrintf ("PADDED_ARRAY_VALUES set_padded_array buffer i_f %d %p %d %d %d: ",
                    i_f,(void*)&padded_face_[ia0],0,iy,iz);
          for (int ix=0; ix<nx; ix++) {
            int i = ia0+ix+ nx*(iy+ ny*iz);
            CkPrintf (" %6.3g",padded_face_[i]);
          }
          CkPrintf ("\n");
        }
      }
    }
#endif      

#ifdef DEBUG_PADDED_ARRAY    
    CkPrintf ("TRACE_PADDED_ARRAY DataMsg::set_padded_array i_f %d %p sum %g\n",
              i_f,&padded_face_[ia0],
              cello::sum(&padded_field[ia0],nx,ny,nz,0,0,0,nx,ny,nz));
#endif    
    // move to next field section in padded array
    //    padded_field += na;
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
  const int na = (n3_pf_[0]*n3_pf_[1]*n3_pf_[2]);
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE na %s:%d %d\n",__FILE__,__LINE__,na);
#endif  
  if (na > 0) {
    size += sizeof(int); // padded_face_field_list_.size()

    int nf = padded_face_field_list_.size();
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE nf %s:%d %d\n",__FILE__,__LINE__,nf);
#endif    
    size += nf*sizeof(int); // padded_face_field_list_
    size += (nf*na)*sizeof(cello_float); // padded_face_

    size += 3*sizeof(int); // ma3_pf_[3]
    size += 3*sizeof(int); // n3_pf_[3]
    size += 3*sizeof(int); // iam3_pf_[3]
    size += 3*sizeof(int); // if3_pf_[3]
  }
  size += 8*sizeof(char); // tag_
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("DEBUG %s DataMsg::data_size %p ma3 %d %d %d\n",
            tag_,(void *)this,ma3_pf_[0],ma3_pf_[1],ma3_pf_[2]);
  CkPrintf ("DEBUG %s DataMsg::data_size %p n3 %d %d %d\n",
            tag_,(void *)this,n3_pf_[0],n3_pf_[1],n3_pf_[2]);
  CkPrintf ("DEBUG %s DataMsg::data_size %p iam3 %d %d %d\n",
            tag_,(void *)this,iam3_pf_[0],iam3_pf_[1],iam3_pf_[2]);
  CkPrintf ("DEBUG %s DataMsg::data_size %p if3 %d %d %d\n",
            tag_,(void *)this,if3_pf_[0],if3_pf_[1],if3_pf_[2]);
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

  const int na = (n3_pf_[0]*n3_pf_[1]*n3_pf_[2]);
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE na %s:%d %d\n",__FILE__,__LINE__,na);
#endif  
  (*pi++) = na;
  if (na > 0) {
    const int nf = padded_face_field_list_.size();
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE nf %s:%d %d\n",__FILE__,__LINE__,nf);
#endif    
    (*pi++) = nf;
    for (int i=0; i<nf; i++) {
      (*pi++) = padded_face_field_list_[i];
    }
    for (int i=0; i<(nf*na); i++) {
      (*pcf++) = padded_face_[i];
    }

    for (int i=0; i<3; i++) {
      (*pi++) = ma3_pf_[i];
      (*pi++) = n3_pf_[i];
      (*pi++) = iam3_pf_[i];
      (*pi++) = if3_pf_[i];
    }
  
#ifdef DEBUG_PADDED_ARRAY
  CkPrintf ("DEBUG %s DataMsg::save_data %p ma3 %d %d %d\n",
            tag_,(void *)this,ma3_pf_[0],ma3_pf_[1],ma3_pf_[2]);
  CkPrintf ("DEBUG %s DataMsg::save_data %p n3 %d %d %d\n",
            tag_,(void *)this,n3_pf_[0],n3_pf_[1],n3_pf_[2]);
  CkPrintf ("DEBUG %s DataMsg::save_data %p iam3 %d %d %d\n",
            tag_,(void *)this,iam3_pf_[0],iam3_pf_[1],iam3_pf_[2]);
  CkPrintf ("DEBUG %s DataMsg::save_data %p if3 %d %d %d\n",
            tag_,(void *)this,if3_pf_[0],if3_pf_[1],if3_pf_[2]);
#endif
  }
  strncpy (pc,tag_,8);
  pc+=8*sizeof(char);
  
#ifdef DEBUG_PADDED_ARRAY  
  print("data_msg:save_data");
#endif  
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

  const int na = (*pi++);
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE na %s:%d %d\n",__FILE__,__LINE__,na);
#endif  
  if (na > 0) {
    const int nf = (*pi++);
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE nf %s:%d %d\n",__FILE__,__LINE__,nf);
#endif    
    padded_face_field_list_.resize(nf);
    for (int i=0; i<nf; i++) {
      padded_face_field_list_[i] = (*pi++);
    }
    padded_face_.resize(nf*na);
    for (int i=0; i<(nf*na); i++) {
      padded_face_[i] = (*pcf++);
    }

    for (int i=0; i<3; i++) {
      ma3_pf_[i]  = (*pi++);
      n3_pf_[i]   = (*pi++);
      iam3_pf_[i] = (*pi++);
      if3_pf_[i]  = (*pi++);
    }

#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG DataMsg::load_data %p ma3 %d %d %d\n",
              (void *)this,ma3_pf_[0],ma3_pf_[1],ma3_pf_[2]);
    CkPrintf ("DEBUG DataMsg::load_data %p n3 %d %d %d\n",
              (void *)this,n3_pf_[0],n3_pf_[1],n3_pf_[2]);
    CkPrintf ("DEBUG DataMsg::load_data %p iam3 %d %d %d\n",
              (void *)this,iam3_pf_[0],iam3_pf_[1],iam3_pf_[2]);
    CkPrintf ("DEBUG DataMsg::load_data %p if3 %d %d %d\n",
              (void *)this,if3_pf_[0],if3_pf_[1],if3_pf_[2]);
#endif
    ASSERT2 ("DataMsg::load_data()",
             "Expecting buffer size %d actual size %d",
             data_size(),(pc-buffer),
             (data_size() == (pc-buffer)));
  }
  strncpy (tag_,pc,8);
  tag_[8] = 0;
  pc+=8*sizeof(char);
    
#ifdef DEBUG_PADDED_ARRAY  
  print("data_msg:load_data");
#endif  
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

  const int na = n3_pf_[0]*n3_pf_[1]*n3_pf_[2];
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE na %s:%d %d\n",__FILE__,__LINE__,na);
#endif  
  if (na>0) {
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG %s DataMsg::update %p ma3 %d %d %d\n",
              tag_,(void *)this,ma3_pf_[0],ma3_pf_[1],ma3_pf_[2]);
    CkPrintf ("DEBUG %s DataMsg::update %p n3 %d %d %d\n",
              tag_,(void *)this,n3_pf_[0],n3_pf_[1],n3_pf_[2]);
    CkPrintf ("DEBUG %s DataMsg::update %p iam3 %d %d %d\n",
              tag_,(void *)this,iam3_pf_[0],iam3_pf_[1],iam3_pf_[2]);
    CkPrintf ("DEBUG %s DataMsg::update %p if3 %d %d %d\n",
              tag_,(void *)this,if3_pf_[0],if3_pf_[1],if3_pf_[2]);
#endif
    const int mx  ( ma3_pf_[0]),  my( ma3_pf_[1]),  mz( ma3_pf_[2]);
    const int nx  (  n3_pf_[0]),  ny(  n3_pf_[1]),  nz(  n3_pf_[2]);
    const int iaxm(iam3_pf_[0]),iaym(iam3_pf_[1]),iazm(iam3_pf_[2]);
    const int ifx ( if3_pf_[0]), ify( if3_pf_[1]), ifz( if3_pf_[2]);
    const int na(mx*my*mz);
    const int nb(nx*ny*nz);
    
    const int nf = padded_face_field_list_.size();
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG_PROLONG_DATA_MSG TRACE nf %s:%d %d\n",__FILE__,__LINE__,nf);
#endif    

#ifdef DEBUG_PADDED_ARRAY
    CkPrintf ("DEBUG_PROLONG_DATA_MSG data %p padded_face_array_allocate %d %d %d  %d*%d*%d*%d\n",
              (void *)data,ifx,ify,ifz,nf,mx,my,mz);
              
#endif

    cello_float * padded_array =
      data->field_data()->padded_array_allocate (ifx,ify,ifz, nf, mx,my,mz);
#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG_PADDED_ARRAY_ALLOCATE data %p %p\n",data,&padded_array[0]);
#endif    
    
    const int ia_start = iaxm + mx*(iaym + my*iazm);

#ifdef DEBUG_PADDED_ARRAY  
    CkPrintf ("DEBUG_PROLONG_DATA_MSG padded_face_.size() %d [%d  %d %d %d]\n",
              padded_face_.size(),nf,nx,ny,nz);
    CkPrintf ("DEBUG_PROLONG_DATA_MSG padded_array %p [%d  %d %d %d]\n",
              padded_array,nf,mx,my,mz);
#endif    
    for (int i_f=0; i_f<nf; i_f++) {
#ifdef DEBUG_PADDED_ARRAY  
      CkPrintf ("DEBUG_PROLONG_DATA_MSG padded_field_list[%d] = %d\n",
                i_f,padded_field_list_[i_f]);
#endif      
      double sum_a=0.0,sum_b=0.0;
      int count=0;
      const int ia0 (na*i_f);
#ifdef TRACE_PADDED_ARRAY_VALUES
      if (i_f == 0) {
        CkPrintf ("PADDED_ARRAY_VALUES %s:%d tag %s i_f %d update array %p padded_array before\n",
                  __FILE__,__LINE__,tag_,i_f,(void*)&padded_array[ia0]);
        for (int iz=0; iz<mz; iz++) {
          for (int iy=0; iy<my; iy++) {
            CkPrintf ("PADDED_ARRAY_VALUES update array pre tag %s i_f %d %p %d %d %d: ",
                      tag_,i_f,(void*)&padded_array[ia0],0,iy,iz);
            for (int ix=0; ix<mx; ix++) {
              int i = ia0+ix+ mx*(iy+ my*iz);
              CkPrintf (" %6.3g",padded_array[i]);
            }
            CkPrintf ("\n");
          }
        }
      }
#endif      
      
      const int ib0 (nb*i_f);
      for (int iz=0; iz<nz; iz++) {
        for (int iy=0; iy<ny; iy++) {
          for (int ix=0; ix<nx; ix++) {
            int ia = ia0 + (ix + mx*(iy + my*iz) + ia_start);
            int ib = ib0 + ix + nx*(iy + ny*iz);
            ASSERT6("DataMsg::update",
                    "padded_array invalid access 0 <= %d (%d : %d %d %d) < %ld",
                    ib,i_f,ix,iy,iz,padded_face_.size(),
                    0 <= ib && ib < padded_face_.size());
            padded_array[ia] = padded_face_[ib];
#ifdef DEBUG_PADDED_ARRAY  
            CkPrintf ("DEBUG_PADDED_ARRAY write %d %d %d [%d] = %g\n",
                      ifx,ify,ifz,ia,padded_array[ia]);
#endif            
            sum_a += padded_array[ia];
            sum_b += padded_face_[ib];
            ++count;
            ASSERT8("DataMsg::update",
                    "nan in sum (%d : %d %d %d) / (%d : %d %d %d)",
                    i_f,ix,iy,iz,nf,nx,ny,nz,
                    (sum_a == sum_a));
          }
        }
      }
#ifdef DEBUG_TRACE_SUMS
    {
      int count = 0;
      double sum = 0.0;
      int zero = 0;
      for (int iz=0; iz<nz; iz++) {
        for (int iy=0; iy<ny; iy++) {
          for (int ix=0; ix<nx; ix++) {
            int i=ib0+ix+nx*(iy+ny*iz);
            count++;
            sum += padded_face_[i];
            if (padded_face_[i] == 0) ++zero;
          }
        }
      }
      CkPrintf ("DEBUG_TRACE_SUM recv buffer [%d %d %d] n3 %d %d %d i_f %d %s %d %d %8.4g\n",
                if3_pf_[0],if3_pf_[1],if3_pf_[2],nx,ny,nz,
                i_f,tag_,zero,count,sum);
    }
    {
      int count = 0;
      double sum = 0.0;
      int zero = 0;
      for (int iz=0; iz<mz; iz++) {
        for (int iy=0; iy<my; iy++) {
          for (int ix=0; ix<mx; ix++) {
            int i=ia0+ix+mx*(iy+my*iz);
            count++;
            sum += padded_array[i];
            if (padded_array[i] == 0) ++zero;
          }
        }
      }
      CkPrintf ("DEBUG_TRACE_SUM recv array extra [%d %d %d] %d n3 %d %d %d i_f %d A %p F %p %d %d %8.4g\n",
                if3_pf_[0],if3_pf_[1],if3_pf_[2],ia0,mx,my,mz,i_f,
                &padded_array[ia0],data->field_data(),zero,count,sum);
    }
#endif    
#ifdef TRACE_PADDED_ARRAY_VALUES
      if (i_f == 0) {
        CkPrintf ("PADDED_ARRAY_VALUES %s:%d tag %s i_f %d update %p padded_face_\n",
                  __FILE__,__LINE__,tag_,i_f,&padded_face_[ib0]);
        for (int iz=0; iz<nz; iz++) {
          for (int iy=0; iy<ny; iy++) {
            CkPrintf ("PADDED_ARRAY_VALUES update buffer post tag %s i_f %d %p %d %d %d: ",
                      tag_,i_f,&padded_face_[ib0],0,iy,iz);
            for (int ix=0; ix<nx; ix++) {
              int i = ib0 + ix+ nx*(iy+ ny*iz);
              CkPrintf (" %6.3g",padded_face_[i]);
            }
            CkPrintf ("\n");
          }
        }
        CkPrintf ("PADDED_ARRAY_VALUES %s:%d update %p padded_array after %d %d %d \n",
                  __FILE__,__LINE__,&padded_array[ia0],mx,my,mz);
        for (int iz=0; iz<mz; iz++) {
          for (int iy=0; iy<my; iy++) {
            CkPrintf ("PADDED_ARRAY_VALUES update array post tag %s i_f %d %p %d %d %d: ",
                      tag_,i_f,&padded_array[ia0],0,iy,iz);
            for (int ix=0; ix<mx; ix++) {
              int i = ia0 + ix+ mx*(iy+ my*iz);
              CkPrintf (" %6.3g",padded_array[i]);
            }
            CkPrintf ("\n");
          }
        }
      }
#endif
#ifdef DEBUG_PADDED_ARRAY    
      CkPrintf ("TRACE_PADDED_ARRAY DataMsg::update i_f %d ia0 %d %p sum %g (total %g)\n",
                i_f,ia0,&padded_array[ia0],
                cello::sum(&padded_array[ia0],mx,my,mz,iaxm,iaym,iazm,nx,ny,nz),
                cello::sum(&padded_array[ia0],mx,my,mz,0,0,0,mx,my,mz));
#endif    
#ifdef DEBUG_PADDED_ARRAY  
      CkPrintf ("DEBUG_PROLONG_DATA_MSG field %d/%d sum array %g sum face %g count %d\n",
                i_f,nf,sum_a,sum_b,count);
#endif      
    }
  }
#ifdef DEBUG_PADDED_ARRAY  
  print("data_msg:update");
#endif  

}

//----------------------------------------------------------------------

void DataMsg::print (const char * message) const
{
#ifdef DEBUG_PADDED_ARRAY  
  CkPrintf ("%s DATA_MSG %p\n",  message,(void *)this);
  CkPrintf ("%s DATA_MSG field_face_    = %p\n",
            message,(void*)field_face_);
  CkPrintf ("%s DATA_MSG field_data_    = %p\n",
            message,(void*)field_data_);
  CkPrintf ("%s DATA_MSG particle_data_ = %p\n",
            message,(void*)particle_data_);
  CkPrintf ("%s DATA_MSG particle_data_delete_ = %d\n",
            message,particle_data_delete_?1:0);
  CkPrintf ("%s DATA_MSG |face_fluxes_list_| = %lu\n",
            message,face_fluxes_list_.size());
  CkPrintf ("%s DATA_MSG |face_fluxes_delete_| = %lu\n",
            message,face_fluxes_delete_.size());
  CkPrintf ("%s DATA_MSG field_face_delete_ = %d\n",
            message,field_face_delete_?1:0);
  CkPrintf ("%s DATA_MSG field_data_delete_ = %d\n",
            message,field_data_delete_?1:0);
  CkPrintf ("%s DATA_MSG padded_face_.sum = %f\n", message,
            std::accumulate(padded_face_.begin(),padded_face_.end(),0.0));
  CkPrintf ("%s DATA_MSG padded_face_field_list_.sum = %d\n", message,
            std::accumulate
            (padded_face_field_list_.begin(),
             padded_face_field_list_.end(),0));
  CkPrintf ("%s DATA_MSG padded ma3_pf_   = %d %d %d\n",
            message,ma3_pf_[0],ma3_pf_[1],ma3_pf_[2]);
  CkPrintf ("%s DATA_MSG padded n3_pf_   = %d %d %d\n",
            message,n3_pf_[0],n3_pf_[1],n3_pf_[2]);
  CkPrintf ("%s DATA_MSG padded iam3_pf_   = %d %d %d\n",
            message,iam3_pf_[0],iam3_pf_[1],iam3_pf_[2]);
  CkPrintf ("%s DATA_MSG padded if3_pf_   = %d %d %d\n",
            message,if3_pf_[0],if3_pf_[1],if3_pf_[2]);
    
  fflush(stdout);
#endif  
}
