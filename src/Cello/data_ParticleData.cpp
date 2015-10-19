// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 11:48:15 PDT 2014
/// @brief    Implementation of the ParticleData class

#include "data.hpp"

//----------------------------------------------------------------------

ParticleData::ParticleData()
{
}

//----------------------------------------------------------------------

void ParticleData::pup (PUP::er &p)
{
}

//----------------------------------------------------------------------

char * ParticleData::attribute_array (ParticleDescr * particle_descr,
				      int it,int ib,int ia)
{

  bool in_range = true;
  if ( !(0 <= it && it < particle_descr->num_types()) )
    in_range = false;
  if ( !(0 <= ib && ib < num_batches(it)) )
    in_range = false;
  if ( !(0 <= ia && ia < particle_descr->num_attributes(it)) )
    in_range = false;

  char * array = NULL;
  if ( in_range ) {
    int offset = particle_descr->attribute_offset(it,ia);
    int align =  attribute_align_[it][ib];
    array = &attribute_array_[it][ib][0] + (offset + align);
  }
  return array;
}

//----------------------------------------------------------------------

int ParticleData::num_batches (int it) const
{
  if ( !(0 <= it && it < (int) attribute_array_.size() )) return 0;

  return attribute_array_[it].size();
}

//----------------------------------------------------------------------

int ParticleData::num_particles 
(ParticleDescr * particle_descr, int it, int ib) const
{
  if ( !(0 <= it && it < particle_descr->num_types()) ) return 0;
  if ( !(0 <= ib && ib < num_batches(it)) ) return 0;

  int bytes_array = attribute_array_[it][ib].size() - 15;
  const int mb = particle_descr->attribute_bytes(it);
  int size = bytes_array / mb;

  return size;
}

//----------------------------------------------------------------------

int ParticleData::num_particles 
(ParticleDescr * particle_descr,int it) const
{
  int nb = num_batches(it);
  int np = 0;
  for (int ib=0; ib<nb; ib++) {
    np += num_particles(particle_descr,it,ib);
  }
  return np;
}

//----------------------------------------------------------------------

int ParticleData::num_particles 
(ParticleDescr * particle_descr) const
{
  int nt = particle_descr->num_types();
  int np = 0;
  for (int it=0; it<nt; it++) {
    np += num_particles(particle_descr,it);
  }
  return np;
}

//----------------------------------------------------------------------

int ParticleData::insert_particles 
(ParticleDescr * particle_descr,
 int it, int np)
{
  int ib0 = num_batches(it) - 1;
  int ip0 = num_particles(particle_descr,it,ib0);
  const int mb = particle_descr->batch_size();
  if (ib0 < 0 || ip0 == mb) {
    ib0++;
    ip0 = 0;
  }
  int ib = ib0;
  int ip = ip0;
  while (np > 0) {
    int np1 = std::min(mb,np);
    if (ip > 0) {
      np1 = std::min(np1,mb-ip);
    }
    if ( ib > num_batches(it) - 1) {
      attribute_array_[it].resize(ib+1);
      attribute_align_[it].resize(ib+1);
    }
    // allocate 15 extra bytes so there is enough room when aligning
    resize_array_(particle_descr,it,ib,ip+np1);
    np -= np1;
    ib++;
    ip=0;
  }

  // return global index of first particle
  return (ib0*mb+ip0);
}

//----------------------------------------------------------------------

void ParticleData::delete_particles 
(ParticleDescr * particle_descr,
 int it, int ib, const bool * mask)
{
  const int na = particle_descr->attribute_bytes(it);
  const int np = num_particles(particle_descr,it,ib);
  int nd=0;
  char * array = attribute_array(particle_descr,it,ib,0);
  for (int ip=0; ip<np; ip++) {
    if (mask[ip]) {
      nd++;
    } else if (nd>0) {
      for (int ia=0; ia<na; ia++) {
	array [(ip-nd)*na+ia] = array [ip*na+ia];
      }
    }
  }
  printf ("nd = %d\n",nd);
  if (nd>0) {
    resize_array_(particle_descr,it,ib,np-nd);
  }
}

//----------------------------------------------------------------------

void ParticleData::split_particles 
(ParticleDescr * particle_descr, 
 int it, int ib, const bool *m,
 ParticleData * particle_data_split)
{
}

//----------------------------------------------------------------------

void ParticleData::compress 
(ParticleDescr * particle_descr,
 int it, int ib, const bool * m)
{
}

//----------------------------------------------------------------------


void ParticleData::resize_array_(ParticleDescr * particle_descr,
				 int it, int ib, int np)
{
  const int mp = particle_descr->attribute_bytes(it);
  
  attribute_array_[it][ib].resize(mp*(np)+15);
  char * array = &attribute_array_[it][ib][0];
  uintptr_t iarray = (uintptr_t) array;
  int defect = (iarray % 16);
  attribute_align_[it][ib] = (defect == 0) ? 0 : 16-defect;
}
