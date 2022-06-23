// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoParticleData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoParticleData class

#include "io.hpp"

//----------------------------------------------------------------------

IoParticleData::IoParticleData() throw ()
  : Io(),
    particle_data_(0),
    particle_index_(0)

{
  // meta_name_.push_back("particle_num_types");
  // meta_name_.push_back("array_size");
  // meta_name_.push_back("num_particles");
  // meta_name_.push_back("ghosts_allocated");
}

//----------------------------------------------------------------------

void IoParticleData::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  Io::pup(p);

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
