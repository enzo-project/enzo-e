// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 11:48:15 PDT 2014
/// @brief    Implementation of the ParticleData class

#include "data.hpp"

//----------------------------------------------------------------------

ParticleData::ParticleData()
  : attribute_array_(),
    attribute_align_(),
    particle_count_()
{
}

//----------------------------------------------------------------------

void ParticleData::pup (PUP::er &p)
{
  p | attribute_array_;
  p | attribute_align_;
  p | particle_count_;
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

  return particle_count_[it][ib];
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

  check_arrays_(particle_descr,__FILE__,__LINE__);

  // find indices of last batch and particle in batch for return value
  int ib_last = num_batches(it) - 1;
  int ip_last = num_particles(particle_descr,it,ib_last);

  int np_left = np;

  const int np_batch = particle_descr->batch_size();

  if (ib_last < 0 || ip_last == np_batch) {
    ib_last++;
    ip_last = 0;
  }
  int ib_this = ib_last;
  int ip_start = ip_last;
  while (np_left > 0) {

    // number of particles to add in this batch
    int np_this = std::min(np_batch,np_left) - ip_start;

    // resize arrays for a new batch if needed
    if ( ib_this > num_batches(it) - 1) {
      attribute_array_[it].resize(ib_this+1);
      attribute_align_[it].resize(ib_this+1);
      particle_count_[it].resize(ib_this+1);
    }

    // allocate particles
    
    resize_array_(particle_descr,it,ib_this,ip_start+np_this);
    particle_count_[it][ib_this] = ip_start+np_this;

    // prepare for next batch

    np_left -= np_this;
    ib_this++;
    ip_start=0;
  }

  // return global index of first particle
  return (ip_last + np_batch*ib_last);
}

//----------------------------------------------------------------------

void ParticleData::delete_particles 
(ParticleDescr * particle_descr,
 int it, int ib, const bool * mask)
{
  check_arrays_(particle_descr,__FILE__,__LINE__);

  const bool interleaved = particle_descr->interleaved(it);

  int npd=0;

  const int na = particle_descr->num_attributes(it);

  const int np = num_particles(particle_descr,it,ib);

  if (interleaved) {

    const int nb = particle_descr->particle_bytes(it);

    for (int ip=0; ip<np; ip++) {
      if (mask[ip]) {
	npd++;
      } else if (npd>0) {
	for (int ia=0; ia<na; ia++) {
	  char * a = attribute_array(particle_descr,it,ib,ia);
	  for (int i=0; i<nb; i++) {
	    a [i + nb*(ip-npd)] = a [i + nb*ip];
	  }
	}
      }
    }
  } else { // ! interleaved

    for (int ip=0; ip<np; ip++) {
      if (mask[ip]) {
	npd++;
      } else if (npd>0) {
	for (int ia=0; ia<na; ia++) {
	  const int nb = particle_descr->attribute_bytes(it,ia);
	  char * a = attribute_array(particle_descr,it,ib,ia);
	  for (int i=0; i<nb; i++) {
	    a [i + nb*(ip-npd)] = a [i + nb*ip];
	  }
	}
      }
    }
  }

  if (npd>0 && interleaved) {
    resize_array_(particle_descr,it,ib,np-npd);
  }
  particle_count_[it][ib] = np-npd;
}

//----------------------------------------------------------------------

void ParticleData::split_particles 
(ParticleDescr * particle_descr, 
 int it, int ib, const bool *mask,
 ParticleData * particle_data_dest)
{
  check_arrays_(particle_descr,__FILE__,__LINE__);
  const int na = particle_descr->num_attributes(it);
  const int np = num_particles(particle_descr,it,ib);

  // Return if no particles deleted
  int nd=0;
  for (int ip=0; ip<np; ip++) {
    if (mask[ip]) nd++;
  }
  if (nd==0) return;

  // insert nd particles

  int i = particle_data_dest->insert_particles(particle_descr,it,nd);

  // copy particles

  const bool interleaved = particle_descr->interleaved(it);
  if (interleaved) {
    for (int ip=0; ip<np; ip++) {
      const int mp = particle_descr->particle_bytes(it);
      if (mask[ip]) {
	int ibc,ipc;
	particle_descr->index(i,&ibc,&ipc);
	for (int ia=0; ia<na; ia++) {
	  char * array = this->attribute_array(particle_descr,it,ib,ia);
	  char * array_dest = 
	    particle_data_dest->attribute_array(particle_descr,it,ibc,ia);
	  int ma = particle_descr->attribute_bytes(it,ia);
	  for (int ib=0; ib<ma; ib++) {
	    array_dest[ipc*ma+ib] = array[ip*ma+ib];
	  }
	}
      }
      i++;
    }
  }  else {
  }

  // delete particles
  delete_particles(particle_descr,it,ib,mask);
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
  const int mp = particle_descr->particle_bytes(it);

  if (!particle_descr->interleaved(it)) {
    np = particle_descr->batch_size();
  }
  attribute_array_[it][ib].resize(mp*(np) + (PARTICLE_ALIGN - 1) );
  char * array = &attribute_array_[it][ib][0];
  uintptr_t iarray = (uintptr_t) array;
  int defect = (iarray % PARTICLE_ALIGN);
  attribute_align_[it][ib] = (defect == 0) ? 0 : PARTICLE_ALIGN-defect;
}

//----------------------------------------------------------------------

void ParticleData::check_arrays_ (ParticleDescr * particle_descr,
		    std::string file, int line) const
{
  size_t nt = particle_descr->num_types();
  ASSERT4 ("ParticleData::check_arrays_",
	  "%s:%d attribute_array_ is size %d < %d",
	   file.c_str(),line,
	   attribute_array_.size(),nt,
	   attribute_array_.size()>=nt);
  ASSERT4 ("ParticleData::check_arrays_",
	  "%s:%d particle_count_ is size %d < %d",
	   file.c_str(),line,
	   particle_count_.size(),nt,
	   particle_count_.size()>=nt);
  ASSERT4 ("ParticleData::check_arrays_",
	  "%s:%d attribute_align_ is size %d < %d",
	   file.c_str(),line,
	   attribute_align_.size(),nt,
	   attribute_align_.size()>=nt);

  for (size_t it=0; it<nt; it++) {
    size_t nb = num_batches(it);

    ASSERT5 ("ParticleData::check_arrays_",
	     "%s:%d attribute_array_[%d] is size %d < %d",
	     file.c_str(),line,
	     it,attribute_array_[it].size(),nb,
	     attribute_array_[it].size()>=nb);
    ASSERT5 ("ParticleData::check_arrays_",
	     "%s:%d particle_count_[%d] is size %d < %d",
	     file.c_str(),line,
	     it,particle_count_[it].size(),nb,
	     particle_count_[it].size()>=nb);
    ASSERT5 ("ParticleData::check_arrays_",
	     "%s:%d attribute_align_[%d] is size %d < %d",
	     file.c_str(),line,
	     it,attribute_align_[it].size(),nb,
	     attribute_align_[it].size()>=nb);
  }
}
