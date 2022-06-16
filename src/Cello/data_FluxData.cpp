// See LICENSE_CELLO file for license and copyright information

/// @file     data_FluxData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-30
/// @brief    Implementation of the FluxData class

#include "data.hpp"

// #define DEBUG_REFRESH

//======================================================================

void FluxData::allocate
( int nx, int ny, int nz,
  std::vector<int> field_list,
  bool single_array,
  std::vector<int> * cx_list,
  std::vector<int> * cy_list,
  std::vector<int> * cz_list)
{
  field_list_ = field_list;
  unsigned nf = field_list.size();
  block_fluxes_.resize(6*nf,nullptr);
  neighbor_fluxes_.resize(6*nf,nullptr);

  const int ix32[3][2] = { {-1, +1},  {0,  0}, {0,  0} };
  const int iy32[3][2] = { { 0,  0}, {-1, +1}, {0,  0} };
  const int iz32[3][2] = { { 0,  0},  {0, 0}, {-1, +1} };

  // Create face flux objects

  for (unsigned i_f=0; i_f<nf; i_f++) {

    const int index_field = field_list[i_f];
    const int cx = cx_list ? (*cx_list)[i_f] : 0;
    const int cy = cy_list ? (*cy_list)[i_f] : 0;
    const int cz = cz_list ? (*cz_list)[i_f] : 0;

    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
        int ix = ix32[axis][face];
        int iy = iy32[axis][face];
        int iz = iz32[axis][face];
        const int i = index_(axis,face,i_f);
        block_fluxes_[i] = new FaceFluxes
          (Face(ix,iy,iz,axis,face),index_field,nx,ny,nz,cx,cy,cz);
        neighbor_fluxes_[i] = new FaceFluxes
          (Face(ix,iy,iz,axis,face),index_field,nx,ny,nz,cx,cy,cz);
      }
    }
  }

  // Allocate flux storage

  if (single_array) {

    // Allocate face fluxes as a single array
    
    size_t size = 0;
    // Determine array size
    for (unsigned i_f=0; i_f<nf; i_f++) {
      for (int axis=0; axis<3; axis++) {
        for (int face=0; face<2; face++) {
          const int i = index_(axis,face,i_f);
          const int m = block_fluxes_[i]->get_size();
          size += 2*m;
        }
      }
    }

    // Allocate storage
    flux_vector_.resize(size);

    // Initialize FaceFlux pointers
    size = 0;
    for (unsigned i_f=0; i_f<nf; i_f++) {
      for (int axis=0; axis<3; axis++) {
        for (int face=0; face<2; face++) {
          const int i = index_(axis,face,i_f);
          const int m = block_fluxes_[i]->get_size();
          block_fluxes_[i] ->set_storage(&flux_vector_.data()[size]);
          size += m;
          neighbor_fluxes_[i]->set_storage(&flux_vector_.data()[size]);
          size += m;
        }
      }
    }
    
  } else {

    // allocate face fluxes individually

    for (unsigned i_f=0; i_f<nf; i_f++) {
      for (int axis=0; axis<3; axis++) {
        for (int face=0; face<2; face++) {
          int i = index_(axis,face,i_f);
          block_fluxes_[i]->allocate_storage();
          neighbor_fluxes_[i]->allocate_storage();
        }
      }
    }

  }
        
}

//----------------------------------------------------------------------

void FluxData::deallocate()
{
  for (unsigned i=0; i<block_fluxes_.size(); i++) {
    delete block_fluxes_[i];
  }
  block_fluxes_.clear();
  for (unsigned i=0; i<neighbor_fluxes_.size(); i++) {
    delete neighbor_fluxes_[i];
  }
  neighbor_fluxes_.clear();
}

//----------------------------------------------------------------------

int FluxData::data_size () const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH %p FluxData::data_size()\n",CkMyPe(),(void*)this);
  fflush(stdout);
#endif  
  int size = 0;

  size += sizeof(int);
  for (unsigned i=0; i<block_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = block_fluxes_[i];
    unsigned n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    size += sizeof(unsigned);
    size += n;
  }

  size += sizeof(unsigned);
  for (unsigned i=0; i<neighbor_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = neighbor_fluxes_[i];
    unsigned n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    size += sizeof(unsigned);
    size += n;
  }

  SIZE_VECTOR_TYPE(size,int,field_list_);
  SIZE_VECTOR_TYPE(size,cello_float,flux_vector_);
  
  return size;
}

//----------------------------------------------------------------------

char * FluxData::save_data (char * buffer) const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH %p FluxData::save_data()\n",CkMyPe(),(void *)this);
  fflush(stdout);
#endif  
  union {
    unsigned *pu;
    char *    pc;
  };

  pc = (char *) buffer;

  (*pu++) = block_fluxes_.size();

  for (unsigned i=0; i<block_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = block_fluxes_[i];
    unsigned n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    (*pu++) = n;
    pc = face_fluxes->save_data(pc);
  }

  (*pu++) = neighbor_fluxes_.size();

  for (unsigned i=0; i<neighbor_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = neighbor_fluxes_[i];
    unsigned n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    (*pu++) = n;
    pc = face_fluxes->save_data(pc);
  }

  SAVE_VECTOR_TYPE(pc,int,field_list_);
  SAVE_VECTOR_TYPE(pc,cello_float,flux_vector_);

  ASSERT2("FluxData::save_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------

char * FluxData::load_data (char * buffer)
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH %p FluxData::load_data()\n",CkMyPe(),(void *)this);
  fflush(stdout);
#endif  
  union {
    unsigned * pu;
    char *     pc;
  };

  pc = (char *) buffer;

  // block_fluxes_

  block_fluxes_.resize(*pu++);

  for (unsigned i=0; i<block_fluxes_.size(); i++) {
    int n = (*pu++);
    FaceFluxes * face_fluxes = (n==0) ? nullptr : new FaceFluxes;
    block_fluxes_[i] = face_fluxes;
    if (face_fluxes != nullptr) {
      pc = face_fluxes->load_data(pc);
    }
  }

  // neighbor_fluxes_

  neighbor_fluxes_.resize(*pu++);

  for (unsigned i=0; i<neighbor_fluxes_.size(); i++) {
    int n = (*pu++);
    FaceFluxes * face_fluxes = (n==0) ? nullptr : new FaceFluxes;
    neighbor_fluxes_[i] = face_fluxes;
    if (face_fluxes != nullptr) {
      pc = face_fluxes->load_data(pc);
    }
  }

  LOAD_VECTOR_TYPE(pc,int,field_list_);
  LOAD_VECTOR_TYPE(pc,cello_float,flux_vector_);

  if (flux_vector_.size() > 0) {

    // Initialize FaceFlux pointers if needed
    const int nf = field_list_.size();
    size_t size = 0;
    for (int i_f=0; i_f<nf; i_f++) {
      for (int axis=0; axis<3; axis++) {
        for (int face=0; face<2; face++) {
          cello_float * flux_data = flux_vector_.data();
          const int i = index_(axis,face,i_f);
          const int flux_size = block_fluxes_[i]->get_size();
          block_fluxes_[i]->set_storage(&flux_data[size]);
          size += flux_size;
          neighbor_fluxes_[i]->set_storage(&flux_data[size]);
          size += flux_size;
        }
      }
    }
  }
  ASSERT2("FluxData::load_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------
  
