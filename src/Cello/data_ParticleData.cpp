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
    array = &attribute_array_[it][ib][0];
    for (int i=0; i<ia; i++) {
      array += particle_descr->attribute_bytes(it,i);
    }
  }
  return array;
}

//----------------------------------------------------------------------

int ParticleData::num_batches (int it) const
{
  if ( !(0 <= it && it < attribute_array_.size() )) return 0;

  return attribute_array_[it].size();
}

//----------------------------------------------------------------------

int ParticleData::num_particles (ParticleDescr * particle_descr, int it, int ib) const
{
  if ( !(0 <= it && it < particle_descr->num_types()) ) return 0;
  if ( !(0 <= ib && ib < num_batches(it)) ) return 0;

  return attribute_array_[it][ib].size();
}

//----------------------------------------------------------------------

int ParticleData::num_particles (ParticleDescr * particle_descr,int it) const
{
  int nb = num_batches(it);
  int np = 0;
  for (int ib=0; ib<nb; ib++) np += num_particles(particle_descr,it,ib);
  return np;
}

//----------------------------------------------------------------------

int ParticleData::num_particles (ParticleDescr * particle_descr) const
{
  int nt = particle_descr->num_types();
  int np = 0;
  for (int it=0; it<nt; it++) {
    np += num_particles(particle_descr);
  }
  return np;
}

//----------------------------------------------------------------------
