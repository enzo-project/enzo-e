// See LICENSE_CELLO file for license and copyright information

/// @file     data_FluxData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-30
/// @brief    Implementation of the FluxData class

#include "data.hpp"

// #define DEBUG_REFRESH

//======================================================================

int FluxData::data_size () const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FluxData::data_size()\n",CkMyPe());
#endif  
  int size = 0;

  size += sizeof(int);
  for (int i=0; i<face_fluxes_block_.size(); i++) {
    FaceFluxes * face_fluxes = face_fluxes_block_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    size += sizeof(int);
    size += n;
  }

  size += sizeof(int);
  for (int i=0; i<face_fluxes_neighbor_.size(); i++) {
    FaceFluxes * face_fluxes = face_fluxes_neighbor_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    size += sizeof(int);
    size += n;
  }
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

  (*pi++) = face_fluxes_block_.size();

  for (int i=0; i<face_fluxes_block_.size(); i++) {
    FaceFluxes * face_fluxes = face_fluxes_block_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    (*pi++) = n;
    pc = face_fluxes->save_data(pc);
  }

  for (int i=0; i<face_fluxes_neighbor_.size(); i++) {
    FaceFluxes * face_fluxes = face_fluxes_neighbor_[i];
    int n = (face_fluxes==nullptr) ? 0 : face_fluxes->data_size();
    (*pi++) = n;
    pc = face_fluxes->save_data(pc);
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

  face_fluxes_block_.resize(*pi++);

  for (int i=0; i<face_fluxes_block_.size(); i++) {
    int n = (*pi++);
    FaceFluxes * face_fluxes = (n==0) ? nullptr : new FaceFluxes;
    face_fluxes_block_[i] = face_fluxes;
    if (face_fluxes != nullptr) {
      pc = face_fluxes->load_data(pc);
    }
  }

  for (int i=0; i<face_fluxes_neighbor_.size(); i++) {
    int n = (*pi++);
    FaceFluxes * face_fluxes = (n==0) ? nullptr : new FaceFluxes;
    face_fluxes_neighbor_[i] = face_fluxes;
    if (face_fluxes != nullptr) {
      pc = face_fluxes->load_data(pc);
    }
  }

  ASSERT2("FluxData::load_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------
  
