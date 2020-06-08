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
  int level, double dt,
  std::vector<int> field_list,
  std::vector<int> * cx_list,
  std::vector<int> * cy_list,
  std::vector<int> * cz_list)
{
  field_list_ = field_list;
  int n_f = field_list.size();
  block_fluxes_.resize(6*n_f,nullptr);
  neighbor_fluxes_.resize(6*n_f,nullptr);

  const int ixm = -1;
  const int ixp = +1;
  const int iym = (ny == 1) ? 0 : -1;
  const int iyp = (ny == 1) ? 0 : +1;
  const int izm = (nz == 1) ? 0 : -1;
  const int izp = (nz == 1) ? 0 : +1;

  const int rx32[3][2] = { {-1, +1},  {0,  0}, {0,  0} };
  const int ry32[3][2] = { { 0,  0}, {-1, +1}, {0,  0} };
  const int rz32[3][2] = { { 0,  0},  {0, 0}, {-1, +1} };
  
  for (int i_f=0; i_f<n_f; i_f++) {
        
    const int index_field = field_list[i_f];

    field_index_map_[index_field] = i_f;

    const int cx = cx_list ? (*cx_list)[i_f] : 0;
    const int cy = cy_list ? (*cy_list)[i_f] : 0;
    const int cz = cz_list ? (*cz_list)[i_f] : 0;

    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
        int rx = rx32[axis][face];
        int ry = ry32[axis][face];
        int rz = rz32[axis][face];
        FaceFluxes * ff;
        ff = new FaceFluxes
          (Face(rx,ry,rz,rx,ry,rz),index_field,nx,ny,nz,level,dt,cx,cy,cz);
        ff->allocate();
        block_fluxes_[index_(axis,face,i_f)] = ff;
        ff = new FaceFluxes
          (Face(rx,ry,rz,rx,ry,rz),index_field,nx,ny,nz,level,dt,cx,cy,cz);
        ff->allocate();
        neighbor_fluxes_[index_(axis,face,i_f)] = ff;
      }
    }
  }
}

//----------------------------------------------------------------------

int FluxData::data_size () const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FluxData::data_size()\n",CkMyPe());
#endif  
  int size = 0;

  size += sizeof(int);
  for (int i=0; i<block_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = block_fluxes_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    size += sizeof(int);
    size += n;
  }

  size += sizeof(int);
  for (int i=0; i<neighbor_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = neighbor_fluxes_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    size += sizeof(int);
    size += n;
  }

  // field_list_ (field map constructed from list)

  size += sizeof(int);
  const int n = (field_list_.size());
  size += n*sizeof(int);
  
  return size;
}

//----------------------------------------------------------------------

char * FluxData::save_data (char * buffer) const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FluxData::save_data()\n",CkMyPe());
#endif  
  union {
    int  * pi;
    char * pc;
  };

  pc = (char *) buffer;

  (*pi++) = block_fluxes_.size();

  for (int i=0; i<block_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = block_fluxes_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    (*pi++) = n;
    pc = face_fluxes->save_data(pc);
  }

  (*pi++) = neighbor_fluxes_.size();

  for (int i=0; i<neighbor_fluxes_.size(); i++) {
    FaceFluxes * face_fluxes = neighbor_fluxes_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    (*pi++) = n;
    pc = face_fluxes->save_data(pc);
  }

  (*pi++) = field_list_.size();
  
  for (int i=0; i<field_list_.size(); i++) {
    (*pi++) = field_list_[i];
  }
  
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
  CkPrintf ("%d DEBUG_REFRESH FluxData::load_data()\n",CkMyPe());
#endif  
  union {
    int  * pi;
    char * pc;
  };

  pc = (char *) buffer;

  // block_fluxes_

  block_fluxes_.resize(*pi++);

  for (int i=0; i<block_fluxes_.size(); i++) {
    int n = (*pi++);
    FaceFluxes * face_fluxes = (n==0) ? nullptr : new FaceFluxes;
    block_fluxes_[i] = face_fluxes;
    if (face_fluxes != nullptr) {
      pc = face_fluxes->load_data(pc);
    }
  }

  // neighbor_fluxes_

  neighbor_fluxes_.resize(*pi++);

  for (int i=0; i<neighbor_fluxes_.size(); i++) {
    int n = (*pi++);
    FaceFluxes * face_fluxes = (n==0) ? nullptr : new FaceFluxes;
    neighbor_fluxes_[i] = face_fluxes;
    if (face_fluxes != nullptr) {
      pc = face_fluxes->load_data(pc);
    }
  }

  // field_list_ (and reconstruct field_index_map_)

  int n = (*pi++);
  field_list_.resize(n);
  for (int i=0; i<field_list_.size(); i++) {
    field_list_[i] = (*pi++);
    field_index_map_[field_list_[i]] = i;
  }

  ASSERT2("FluxData::load_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------
  
