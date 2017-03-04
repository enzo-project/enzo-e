// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 11:48:15 PDT 2014
/// @brief    Implementation of the ParticleData class

#include "data.hpp"
#include <algorithm>

// #define DEBUG_PARTICLES

int64_t ParticleData::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

ParticleData::ParticleData()
  : attribute_array_(),
    attribute_align_(),
    particle_count_()
{
  ++counter[cello::index_static()]; 
}

//----------------------------------------------------------------------

bool ParticleData::operator== (const ParticleData & particle_data) throw ()
{
  return 
    (attribute_array_ == particle_data.attribute_array_) &&
    (attribute_align_ == particle_data.attribute_align_) &&
    (particle_count_  == particle_data.particle_count_);
}

//----------------------------------------------------------------------

void ParticleData::pup (PUP::er &p)
{
  p | attribute_array_;
  p | attribute_align_;
  p | particle_count_;
}

//----------------------------------------------------------------------

ParticleData::~ParticleData()
{
  --counter[cello::index_static()]; 
}
//----------------------------------------------------------------------

char * ParticleData::attribute_array (ParticleDescr * particle_descr,
				      int it,int ia,int ib)
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
(ParticleDescr * particle_descr) const
{
  const int nt = particle_descr->num_types();
  int np = 0;
  for (int it=0; it<nt; it++) {
    np += num_particles(particle_descr,it);
  }
  return np;
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
(ParticleDescr * particle_descr, int it, int ib) const
{
  if ( !(0 <= it && it < particle_descr->num_types()) ) return 0;
  if ( !(0 <= ib && ib < num_batches(it)) ) return 0;

  return particle_count_[it][ib];
}

//----------------------------------------------------------------------

int ParticleData::insert_particles 
(ParticleDescr * particle_descr,
 int it, int np)
{
  if (np==0) return 0;

  check_arrays_(particle_descr,__FILE__,__LINE__);

  // find indices of last batch and particle in batch for return value

  const int nb = num_batches(it);
  const int mb = particle_descr->batch_size();

  int ib_last,ip_last;

  if (nb == 0) {
    ib_last = 0;
    ip_last = 0;
  } else {
    ib_last = nb - 1;
    ip_last = num_particles(particle_descr,it,ib_last);
    if (ip_last % mb == 0) {
      ib_last++;
      ip_last = 0;
    }
  }

  int ib_this = ib_last;
  int ip_start = ip_last;

  int np_left = np;

  while (np_left > 0) {

    // number of particles to add in this batch
    int np_this = std::min(mb-ip_start,np_left);

    // resize arrays for a new batch if needed
    if ( ip_start == 0) {

      ASSERT1("ParticleData::insert_particles",
	     "Trying to insert negative particles: ib_this = %d",
	      ib_this, ib_this >= 0);

      attribute_array_[it].resize(ib_this+1);
      attribute_align_[it].resize(ib_this+1);
      particle_count_ [it].resize(ib_this+1);
    }

    // allocate particles
    
    resize_attribute_array_(particle_descr,it,ib_this,ip_start+np_this);

    // prepare for next batch

    np_left -= np_this;
    ib_this++;
    ip_start=0;
  }

  // return global index of first particle
  return (ip_last + mb*ib_last);
}

//----------------------------------------------------------------------

int ParticleData::delete_particles 
(ParticleDescr * particle_descr,
 int it, int ib, const bool * mask)
{
  check_arrays_(particle_descr,__FILE__,__LINE__);

  const bool interleaved = particle_descr->interleaved(it);

  // number of particle deleted so far
  int npd=0;

  const int na = particle_descr->num_attributes(it);
  const int np = num_particles(particle_descr,it,ib);
  int mp = particle_descr->particle_bytes(it);

  // for each particle
  for (int ip=0; ip<np; ip++) {
    // ... move to next in source if deleting
    if ((mask==NULL) || mask[ip]) {
      npd++;
    } else if (npd>0) {
      // ... else copy the particle attributes back to first opening
      for (int ia=0; ia<na; ia++) {
	if (!interleaved) 
	  mp = particle_descr->attribute_bytes(it,ia);
	const int ny = particle_descr->attribute_bytes(it,ia);
	char * a = attribute_array(particle_descr,it,ia,ib);
	for (int iy=0; iy<ny; iy++) {
	  const int i_old = iy + mp*ip;
	  const int i_new = iy + mp*(ip-npd);
	  a [i_new] = a [i_old];
	}
      }
    }
  }

  if (npd>0) {
    resize_attribute_array_(particle_descr,it,ib,np-npd);
  }

  return npd;
}

//----------------------------------------------------------------------

void ParticleData::scatter 
(ParticleDescr * particle_descr,
 int it, int ib,
 int np, const bool * mask, const int * index,
 int n, ParticleData * particle_array[])
{
  // count number of particles in each particle_array element

  int np_array[n] = {0};

  for (int ip=0; ip<np; ip++) {
    if ((mask == NULL) || mask[ip]) {
      ++np_array[index[ip]];
    }
  }
  
  // insert uninitialized particles
  std::map<ParticleData *, int>  i_array;
  std::map<ParticleData *, bool> is_first;
  for (int k=0; k<n; k++) {
    ParticleData * pd = particle_array[k];
    is_first[pd] = true;
    i_array[pd] = 0;
  }

  for (int k=0; k<n; k++) {
    ParticleData * pd = particle_array[k];
    if (np_array[k]>0 && pd) {
      int i0 = pd->insert_particles (particle_descr,it,np_array[k]);
      if (is_first[pd]) i_array[pd] = i0;
      is_first[pd] = false;
    }
  }

  const bool interleaved = particle_descr->interleaved(it);
  const int na = particle_descr->num_attributes(it);
  int mp = particle_descr->particle_bytes(it);
  const int ib_src = ib;

  int count=0;
  for (int ip_src=0; ip_src<np; ip_src++) {

    if ((mask == NULL) || mask[ip_src]) {
      ++count;
      int k = index[ip_src];
      ParticleData * pd = particle_array[k];
      int i_dst = i_array[pd]++;
      int ib_dst,ip_dst;
      particle_descr->index(i_dst,&ib_dst,&ip_dst);
      for (int ia=0; ia<na; ia++) {
	if (!interleaved) 
	  mp = particle_descr->attribute_bytes(it,ia);
	int ny = particle_descr->attribute_bytes(it,ia);
	char * a_src = attribute_array
	  (particle_descr,it,ia,ib_src);
	char * a_dst = pd->attribute_array
	  (particle_descr,it,ia,ib_dst);
	for (int iy=0; iy<ny; iy++) {
	  a_dst [iy + mp*ip_dst] = a_src [iy + mp*ip_src];
	}
      }
    }

  }
}

//----------------------------------------------------------------------

int ParticleData::gather 
(ParticleDescr * particle_descr, int it, 
 int n, ParticleData * particle_array[])
{
  int count = 0;
  
  // Sort particle array to simplify skipping duplicates
  ParticleData * particle_array_sorted[n];
  for (int i=0; i<n; i++) particle_array_sorted[i] = particle_array[i];
  std::sort(&particle_array_sorted[0],
	    &particle_array_sorted[n]);

  // count number of particles to insert
  int np = 0;
  for (int k=0; k<n; k++) {
    // ... skipping duplicate ParticleData objects
    if (k>0 && (particle_array_sorted[k] == 
		particle_array_sorted[k-1])) continue;
    ParticleData * pd = particle_array_sorted[k];
    np += pd ? pd->num_particles(particle_descr,it) : 0;
  }

  // insert uninitialized particles
  int i_dst = insert_particles(particle_descr,it,np);

  int ib_dst,ip_dst;
  particle_descr->index (i_dst,&ib_dst,&ip_dst);

  // initialize particles
  const bool interleaved = particle_descr->interleaved(it);
  const int na = particle_descr->num_attributes(it);
  int mp = particle_descr->particle_bytes(it);

  const int mb = particle_descr->batch_size();

  for (int k=0; k<n; k++) {
    // ...skip duplicate ParticleData objects
    if (k>0 && (particle_array_sorted[k] == 
		particle_array_sorted[k-1])) continue;
    ParticleData * pd = particle_array_sorted[k];
    const int nb = pd ? pd->num_batches(it) : 0;
    for (int ib=0; ib<nb; ib++) {
      const int np = pd->num_particles(particle_descr,it,ib);
      count += np;
      for (int ip=0; ip<np; ip++) {
	for (int ia=0; ia<na; ia++) {
	  if (!interleaved) 
	    mp = particle_descr->attribute_bytes(it,ia);
	  int ny = particle_descr->attribute_bytes(it,ia);
	  char * a_src = pd->attribute_array
	    (particle_descr,it,ia,ib);
	  char * a_dst = attribute_array(particle_descr,it,ia,ib_dst);
	  for (int iy=0; iy<ny; iy++) {
	    a_dst [iy + mp*ip_dst] = a_src [iy + mp*ip];
	  }
	}
	ip_dst = (ip_dst + 1) % mb;
	if (ip_dst == 0) ib_dst++;
      }
    }
  }
  return count;
}

//----------------------------------------------------------------------

void ParticleData::compress (ParticleDescr * particle_descr)
{
  const int nt = particle_descr->num_types();
  for (int it=0; it<nt; it++) {
    compress (particle_descr,it);
  }
}

//----------------------------------------------------------------------

void ParticleData::compress (ParticleDescr * particle_descr, int it)
{
  const int nb = num_batches(it);
  const int mb = particle_descr->batch_size();
  const int na = particle_descr->num_attributes(it);

  const bool interleaved = particle_descr->interleaved(it);

  // find first batch with space in it

  int ib_src, ip_src; // source batch and particle indices
  int ib_dst, ip_dst; // destination batch and particle indices

  int np_dst; // number of particles in ib_dst batch
  int np_src; // number of particles in ib_src batch

  // find destination: first empty spot 

  ib_dst=0;
  np_dst=(ib_dst<nb) ? num_particles(particle_descr,it,ib_dst) : 0;
  while (ib_dst<nb && np_dst == mb) {
    ib_dst++;
  } // assert ib_dst == nb || np_dst < mb
  ip_dst = np_dst;

  // first source: next non-empty spot
  ib_src = ib_dst + 1;
  ip_src = 0;
  np_src = num_particles(particle_descr,it,ib_src);

  resize_attribute_array_ (particle_descr,it,ib_dst,mb);
  np_dst = mb;

  int mp = particle_descr->particle_bytes(it);

  while (ib_src < nb && ip_src < np_src) {

    for (int ia=0; ia<na; ia++) {
      if (!interleaved) {
	mp = particle_descr->attribute_bytes(it,ia);
      }
      int ny = particle_descr->attribute_bytes(it,ia);
      char * a_src = attribute_array(particle_descr,it,ia,ib_src);
      char * a_dst = attribute_array(particle_descr,it,ia,ib_dst);
      for (int iy=0; iy<ny; iy++) {
	a_dst [iy + mp*ip_dst] = a_src [iy + mp*ip_src];
      }
    }

    ip_dst++;
    if (ip_dst>=np_dst) {
      ip_dst = 0;
      ib_dst++;
      np_dst = mb;
      if (ib_dst < nb) resize_attribute_array_ (particle_descr,it,ib_dst,mb);
    }

    ip_src++;
    if (ip_src>=np_src) {
      ip_src = 0;
      ib_src++;
      if (ib_src < nb) np_src = num_particles(particle_descr,it,ib_src);
    }

  }

  // deallocate empty batches?
}
  

//----------------------------------------------------------------------

float ParticleData::efficiency (ParticleDescr * particle_descr)
{
  int64_t bytes_min=0,bytes_used=0;
  const int mb = particle_descr->batch_size();

  const int nt = particle_descr->num_types();
  for (int it=0; it<nt; it++) {
    const int nb = num_batches(it);
    const int mp = particle_descr->particle_bytes(it);
    for (int ib=0; ib<nb; ib++) {
      const int np = num_particles(particle_descr,it,ib);
      bytes_min += np*mp;
      bytes_used += mb*mp;
    }
  }
  return (bytes_used) ? 1.0*bytes_min/bytes_used : 1.0;

}

//----------------------------------------------------------------------

float ParticleData::efficiency (ParticleDescr * particle_descr, int it)
{
  int64_t bytes_min=0,bytes_used=0;
  const int mb = particle_descr->batch_size();

  const int nb = num_batches(it);
  const int mp = particle_descr->particle_bytes(it);
  for (int ib=0; ib<nb; ib++) {
    const int np = num_particles(particle_descr,it,ib);
    bytes_min += np*mp;
    bytes_used += mb*mp;
  }
  return (bytes_used) ? 1.0*bytes_min/bytes_used : 1.0;
}

//----------------------------------------------------------------------

float ParticleData::efficiency (ParticleDescr * particle_descr, int it, int ib)
{

  const int mp = particle_descr->particle_bytes(it);
  const int np = num_particles(particle_descr,it,ib);
  const int mb = particle_descr->batch_size();

  const int64_t bytes_min  = (int64_t)np*mp;
  const int64_t bytes_used = (int64_t)mb*mp;

  return (bytes_used) ? 1.0*bytes_min/bytes_used : -1.0;
  
}

//----------------------------------------------------------------------

bool ParticleData::position 
(
 ParticleDescr * particle_descr,
 int it, int ib,
 double * x, double * y, double * z)
{
  const int ia_x = particle_descr->attribute_position(it,0);
  const int ia_y = particle_descr->attribute_position(it,1);
  const int ia_z = particle_descr->attribute_position(it,2);

  bool l_return = false;
  if (x && (ia_x != -1)) {
    const int type_x = particle_descr->attribute_type(it,ia_x);
    l_return = true;
    if (cello::type_is_float(type_x)) {
      copy_attribute_float_ (particle_descr,type_x,it,ib,ia_x,x);
    } else if (cello::type_is_int(type_x)) {
      copy_position_int_ (particle_descr,type_x,it,ib,ia_x,x);
    }
  }

  if (y && (ia_y != -1)) {
    const int type_y = particle_descr->attribute_type(it,ia_y);
    l_return = true;
    if (cello::type_is_float(type_y)) {
      copy_attribute_float_ (particle_descr,type_y,it,ib,ia_y,y);
    } else if (cello::type_is_int(type_y)) {
      copy_position_int_ (particle_descr,type_y,it,ib,ia_y,y);
    }
  }

  if (z && (ia_z != -1)) {
    const int type_z = particle_descr->attribute_type(it,ia_z);
    l_return = true;
    if (cello::type_is_float(type_z)) {
      copy_attribute_float_ (particle_descr,type_z,it,ib,ia_z,z);
    } else if (cello::type_is_int(type_z)) {
      copy_position_int_ (particle_descr,type_z,it,ib,ia_z,z);
    }
  }

  return l_return;
}

//----------------------------------------------------------------------

bool ParticleData::velocity 
(
 ParticleDescr * particle_descr,
 int it, int ib,
 double * vx, double * vy, double * vz)
{
  // If velocitys are floating-point, return the requested values directly
  int ia_x = particle_descr->attribute_velocity(it,0);
  int ia_y = particle_descr->attribute_velocity(it,1);
  int ia_z = particle_descr->attribute_velocity(it,2);
  bool l_return = false;
  if (vx && (ia_x != -1)) {
    const int type_x = particle_descr->attribute_type(it,ia_x);
    l_return = true;
    copy_attribute_float_ (particle_descr,type_x,it,ib,ia_x,vx);
  }

  if (vy && (ia_y != -1)) {
    const int type_y = particle_descr->attribute_type(it,ia_y);
    l_return = true;
    copy_attribute_float_ (particle_descr,type_y,it,ib,ia_y,vy);
  }

  if (vz && (ia_z != -1)) {
    const int type_z = particle_descr->attribute_type(it,ia_z);
    l_return = true;
    copy_attribute_float_ (particle_descr,type_z,it,ib,ia_z,vz);
  }

  return l_return;
}

//----------------------------------------------------------------------

void ParticleData::position_update 
(ParticleDescr * particle_descr,int it, int ib, 
 long double dx, long double dy, long double dz)
{
  const int ia_x = particle_descr->attribute_position(it,0);
  const int ia_y = particle_descr->attribute_position(it,1);
  const int ia_z = particle_descr->attribute_position(it,2);

  if (dx != 0.0 && (ia_x != -1)) {
    const int type_x = particle_descr->attribute_type(it,ia_x);
    if (cello::type_is_float(type_x)) {
      update_attribute_float_ (particle_descr,type_x,it,ib,ia_x,dx);
    } else if (cello::type_is_int(type_x)) {
      update_position_int_ (particle_descr,type_x,it,ib,ia_x,dx);
    }
  }
  if (dy != 0.0 && (ia_y != -1)) {
    const int type_y = particle_descr->attribute_type(it,ia_y);
    if (cello::type_is_float(type_y)) {
      update_attribute_float_ (particle_descr,type_y,it,ib,ia_y,dy);
    } else if (cello::type_is_int(type_y)) {
      update_position_int_ (particle_descr,type_y,it,ib,ia_y,dy);
    }
  }
  if (dz != 0.0 && (ia_z != -1)) {
    const int type_z = particle_descr->attribute_type(it,ia_z);
    if (cello::type_is_float(type_z)) {
      update_attribute_float_ (particle_descr,type_z,it,ib,ia_z,dz);
    } else if (cello::type_is_int(type_z)) {
      update_position_int_ (particle_descr,type_z,it,ib,ia_z,dz);
    }
  }
  
}

//----------------------------------------------------------------------

void ParticleData::copy_attribute_float_ 
(ParticleDescr * particle_descr,
 int type, int it, int ib, int ia, double * coord)
{

  const int dx = particle_descr->stride(it,ia);
  const char * array = attribute_array(particle_descr,it,ia,ib);
  const int np = num_particles(particle_descr,it,ib);
  if (type == type_float) {
    const float * array_f = (float *) array;
    for (int ip=0; ip<np; ip++) coord[ip] = array_f[ip*dx];
  } else if (type == type_double) {
    const double * array_d = (double *) array;
    for (int ip=0; ip<np; ip++) coord[ip] = array_d[ip*dx];
  } else if (type == type_quadruple) {
    const long double * array_q = (long double *) array;
    for (int ip=0; ip<np; ip++) coord[ip] = array_q[ip*dx];
  } else {
    ERROR1("ParticleData::copy_attribute_float_()",
	   "Unknown particle attribute type %d", type);
  }
}

//----------------------------------------------------------------------

void ParticleData::update_attribute_float_ 
(ParticleDescr * particle_descr,
 int type, int it, int ib, int ia, long double da)
{
  const int dx = particle_descr->stride(it,ia);
  char * array = attribute_array(particle_descr,it,ia,ib);
  const int np = num_particles(particle_descr,it,ib);
  if (type == type_float) {
    float * array_f = (float *) array;
    for (int ip=0; ip<np; ip++) array_f[ip*dx] += da;
  } else if (type == type_double) {
    double * array_d = (double *) array;
    for (int ip=0; ip<np; ip++) array_d[ip*dx] += da;
  } else if (type == type_quadruple) {
    long double * array_q = (long double *) array;
    for (int ip=0; ip<np; ip++) array_q[ip*dx] += da;
  } else {
    ERROR1("ParticleData::copy_attribute_float_()",
	   "Unknown particle attribute type %d",
	   type);
  }
}

//----------------------------------------------------------------------

void ParticleData::copy_position_int_ 
(ParticleDescr * particle_descr,
 int type, int it, int ib, int ia, double * coord)
{
  const int dx = particle_descr->stride(it,ia);
  const char * array = attribute_array(particle_descr,it,ia,ib);
  const int np = num_particles(particle_descr,it,ib);
  if (type == type_int8) {
    const int8_t * array_8 = (int8_t *) array;
    for (int ip=0; ip<np; ip++) coord[ip] = double(1.0)*array_8[ip*dx]/PMAX_8;
  } else if (type == type_int16) {
    const int16_t * array_16 = (int16_t *) array;
    for (int ip=0; ip<np; ip++) coord[ip] = double(1.0)*array_16[ip*dx]/PMAX_16;
  } else if (type == type_int32) {
    const int32_t * array_32 = (int32_t *) array;
    for (int ip=0; ip<np; ip++) coord[ip] = double(1.0)*array_32[ip*dx]/PMAX_32;
  } else if (type == type_int64) {
    const int64_t * array_64 = (int64_t *) array;
    for (int ip=0; ip<np; ip++) coord[ip] = double(1.0)*array_64[ip*dx]/PMAX_64;
  } else {
    ERROR1("ParticleData::copy_attribute_float_()",
	   "Unknown particle attribute type %d",
	   type);
  }

}

//----------------------------------------------------------------------

void ParticleData::update_position_int_ 
(ParticleDescr * particle_descr,
 int type, int it, int ib, int ia, int64_t da)
{
  const int dx = particle_descr->stride(it,ia);
  char * array = attribute_array(particle_descr,it,ia,ib);
  const int np = num_particles(particle_descr,it,ib);
  if (type == type_int8) {
    const int8_t da_8 = (int8_t) da / (PMAX_64/PMAX_8);
    int8_t * array_8 = (int8_t *) array;
    for (int ip=0; ip<np; ip++) array_8[ip*dx] += da_8;
  } else if (type == type_int16) {
    const int16_t da_16 = da / (PMAX_64/PMAX_16);
    int16_t * array_16 = (int16_t *) array;
    for (int ip=0; ip<np; ip++) array_16[ip*dx] += da_16;
  } else if (type == type_int32) {
    const int32_t da_32 = da / (PMAX_64/PMAX_32);
    int32_t * array_32 = (int32_t *) array;
    for (int ip=0; ip<np; ip++) array_32[ip*dx] += da_32;
  } else if (type == type_int64) {
    const int64_t da_64 = da;
    int64_t * array_64 = (int64_t *) array;
    for (int ip=0; ip<np; ip++) array_64[ip*dx] += da_64;
  } else {
    ERROR1("ParticleData::copy_attribute_float_()",
	   "Unknown particle attribute type %d",
	   type);
  }
}

//----------------------------------------------------------------------

int ParticleData::data_size (ParticleDescr * particle_descr) const
{

  int size = 0;
  const int nt = particle_descr->num_types();

  // array lengths

  size += sizeof(int); 

  for (int it=0; it<nt; it++) {

    // array[it] lengths

    size += sizeof(int); 

    const int nb = num_batches(it);

    // attribute_align_[it] values
    size += nb*sizeof(char);

    // particle_count_[it] values
    size += nb*sizeof(int);

    for (int ib=0; ib<nb; ib++) {

      // array[it][ib] length
      size += sizeof(int); 

      // attribute_array_[it][ib] values
      const int mp = attribute_array_[it][ib].size();
      size += mp * sizeof(char);
    }
  }
  return size;
}

//----------------------------------------------------------------------

char * ParticleData::save_data (ParticleDescr * particle_descr,
				char * buffer) const
{
  union {
    int  * pi;
    char * pc;
  };

  // NOTE: integers stored first, then char's, to avoid alignment issues

  pc = (char *) buffer;

  //--------------------
  // Store array sizes
  //--------------------

  // ...store number of types

  const int nt = (*pi++) = particle_descr->num_types();

  for (int it=0; it<nt; it++) {

    // ...store number of batches for the type

    int nb = (*pi++) = num_batches(it);

    for (int ib=0; ib<nb; ib++) {

      // ...store particle attribute array lengths

      (*pi++) = attribute_array_[it][ib].size();
      
    }
  }

  // store int particle_count_[it][ib]

  for (int it=0; it<nt; it++) {
    const int nb = particle_count_[it].size();
    for (int ib=0; ib<nb; ib++) {
      (*pi++) = particle_count_[it][ib];
    }
  }

  // store char attribute_align_[it][ib] array

  for (int it=0; it<nt; it++) {
    const int nb = attribute_align_[it].size();
    for (int ib=0; ib<nb; ib++) {
      (*pc++) = attribute_align_[it][ib];
    }
  }

  // store char attribute_array_[it][ib][k]

  for (int it=0; it<nt; it++) {
    const int nb = attribute_array_[it].size();
    for (int ib=0; ib<nb; ib++) {
      const int n = attribute_array_[it][ib].size();
      for (int k=0; k<n; k++) {
	(*pc++) = attribute_array_[it][ib][k];
      }
    }
  }

  return pc;
}

//----------------------------------------------------------------------

char * ParticleData::load_data (ParticleDescr * particle_descr,
				char * buffer)
{
  // NOTE: integers stored first, then char's, to avoid alignment issues

  union {
    int  * pi;
    char * pc;
  };

  pc = (char *) buffer;

  //-----------------------------------------
  // Load array sizes and pre-allocate arrays
  //-----------------------------------------

  // ...load number of types

  const int nt = (*pi++);

  // ... allocate arrays

  ASSERT1("ParticleData::load_data",
	  "Trying to allocate negative particle types: nt = %d",
	  nt, nt >= 0);
  
  attribute_array_.resize(nt);
  attribute_align_.resize(nt);
  particle_count_.resize(nt);

  for (int it=0; it<nt; it++) {

    // ...load number of batches for the type

    int nb = (*pi++);

    // ... allocate arrays[it]

    ASSERT1("ParticleData::load_data",
	    "Trying to allocate negative particle batches: nb = %d",
	    nb, nb >= 0);

    attribute_array_[it].resize(nb);
    attribute_align_[it].resize(nb);
    particle_count_[it].resize(nb);

    for (int ib=0; ib<nb; ib++) {

      // ...load particle attribute array lengths

      const int np = (*pi++);

      ASSERT1("ParticleData::load_data",
	      "Trying to allocate negative particles: np = %d",
	      np, np >= 0);
      
      attribute_array_[it][ib].resize(np);
      
    }
  }

  // load int particle_count_[it][ib]

  for (int it=0; it<nt; it++) {
    const int nb = particle_count_[it].size();
    for (int ib=0; ib<nb; ib++) {
      particle_count_[it][ib] = (*pi++);
    }
  }

  // load char attribute_align_[it][ib] array

  for (int it=0; it<nt; it++) {
    const int nb = attribute_align_[it].size();
    for (int ib=0; ib<nb; ib++) {
      attribute_align_[it][ib] = (*pc++);
    }
  }

  // load char attribute_array_[it][ib][k]

  for (int it=0; it<nt; it++) {
    const int nb = attribute_array_[it].size();
    for (int ib=0; ib<nb; ib++) {
      const int n = attribute_array_[it][ib].size();
      for (int k=0; k<n; k++) {
	attribute_array_[it][ib][k] = (*pc++);
      }
    }
  }

  return pc;
}

//----------------------------------------------------------------------

void ParticleData::debug (ParticleDescr * particle_descr)
{
  const int nt = particle_descr->num_types();

  for (int it=0; it<nt; it++) {
    int nb = num_batches(it);
    int na = particle_descr->num_attributes(it);
    int np = num_particles(particle_descr,it);
    std::string name = particle_descr->type_name(it);

    printf ("particle %p: type name         %d %s\n",   this,it,name.c_str());
    printf ("particle %p:    num_attributes %d %d\n",   this,it,na);
    for (int ia=0; ia<na; ia++)
      printf ("particle %p:       %s\n",
	      this,particle_descr->attribute_name(it,ia).c_str());
    printf ("particle %p:    num_particles  %d %d\n",   this,it,np);
    printf ("particle %p:    num_batches    %d %d\n",   this,it,nb);
    for (int ib=0; ib<nb; ib++) {
      printf ("particle %p:      ",this);
      for (int ip=0; ip<np; ip++) {
	for (int ia=0; ia<na; ia++) {
	  int d = particle_descr->stride(it,ia);
	  char * a = attribute_array(particle_descr,it,ia,ib);
	  int type = particle_descr->attribute_type (it,ia);
	  if (type==type_double) {
	    printf ("%10.5f ",((double*)a)[ip*d]);
	  } else if (type==type_float) {
	    printf ("%10.5f ",((float*)a)[ip*d]);
	  } else if (type==type_int32) {
	    printf ("%6d ",((int32_t*)a)[ip*d]);
	  } else if (type==type_int64) {
	    printf ("%6ld ",((int64_t*)a)[ip*d]);
	  } else {
	    ERROR1("ParticleData::debug()",
		   "Unknown particle attribute type %d",
		   type);
	  }
	}
	printf ("\n");
      }
    }
    
  }
  
}


//----------------------------------------------------------------------

void ParticleData::write_ifrite (ParticleDescr * particle_descr,
				 int it, std::string file_name,
				 double xm, double ym, double zm,
				 double xp, double yp, double zp)
{
  FILE * fp = fopen (file_name.c_str(),"w");

  fprintf (fp,"%d\n",num_particles(particle_descr,it));
  fprintf (fp,"%f %f %f %f %f %f\n",xm,ym,zm,xp,yp,zp);
  const int nb = num_batches(it);
  const int ia_x = particle_descr->attribute_position(it,0);
  const int d = particle_descr->stride(it,ia_x);
  for (int ib=0; ib<nb; ib++) {
    const int np = num_particles(particle_descr,it,ib);
    double x[np], y[np], z[np];
    position (particle_descr,it,ib,x,y,z);
    for (int ip=0; ip<np; ip++) {
      fprintf (fp,"%f %f %f\n",x[ip*d],y[ip*d],z[ip*d]);
    }
  }
  fflush(fp);
  fclose (fp);
}

//======================================================================


void ParticleData::resize_attribute_array_
(ParticleDescr * particle_descr,int it, int ib, int np)
{
  // store number of particles allocated
  particle_count_[it][ib] = np;

  const int mp = particle_descr->particle_bytes(it);

  if (!particle_descr->interleaved(it)) {
    np = particle_descr->batch_size();
  }

  size_t new_size = mp*(np) + (PARTICLE_ALIGN - 1) ;

  if (attribute_array_[it][ib].size() != new_size) {

    ASSERT1("ParticleData::resize_attribute_array_",
	    "Trying to allocate negative particles: new_size = %d",
	    new_size, new_size >= 0);
      
    attribute_array_[it][ib].resize(new_size);
    char * array = &attribute_array_[it][ib][0];
    uintptr_t iarray = (uintptr_t) array;
    int defect = (iarray % PARTICLE_ALIGN);
    attribute_align_[it][ib] = (defect == 0) ? 0 : PARTICLE_ALIGN-defect;
  }
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
