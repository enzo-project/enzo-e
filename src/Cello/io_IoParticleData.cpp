// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoParticleData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoParticleData class

#include "io.hpp"

//----------------------------------------------------------------------

IoParticleData::IoParticleData(const ParticleDescr * particle_descr) throw ()
  : Io(1),
    particle_descr_(0),
    particle_data_(0),
    particle_index_(0)

{
  meta_name_.push_back("particle_num_types");
  meta_name_.push_back("array_size");
  meta_name_.push_back("num_particles");
  meta_name_.push_back("ghosts_allocated");
}

//----------------------------------------------------------------------

void IoParticleData::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  Io::pup(p);

  //  if (p.isUnpacking()) particle_descr_ = new ParticleDescr;
  //  p | *particle_descr_;
  WARNING ("IoParticleData::pup","skipping particle_descr_");
  
  //  if (p.isUnpacking()) particle_data_ = new ParticleData;
  WARNING ("IoParticleData::pup","skipping particle_data_");
  //  p | *particle_data_;
  p | particle_index_;
}

//----------------------------------------------------------------------

void IoParticleData::meta_value
(int index,
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//----------------------------------------------------------------------

void IoParticleData::field_array
(int index,
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
}

//----------------------------------------------------------------------

void IoParticleData::particle_array
(int it, int ib, int ia,
 void ** buffer, std::string * name, int * type,
 int * n, int * k) throw()
{
  Particle particle (particle_descr_,particle_data_);

  if (buffer) (*buffer) = (void * ) 
   		particle.attribute_array(it,ia,ib);
  if (name) {
    char buffer [80];
    sprintf (buffer,"particle %s %s %d",
  	     particle.type_name(it).c_str(),
	     particle.attribute_name(it,ia).c_str(),
	     ib);
    *name = buffer;
  }

  if (type) (*type) = particle.attribute_type(it,ia);

  if (n) (*n) = particle.num_particles (it,ib);
  if (k) (*k) = particle.stride(it,ia);
}
//----------------------------------------------------------------------
