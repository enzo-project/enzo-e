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
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//----------------------------------------------------------------------

void IoParticleData::data_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd,
 int * nx,  int * ny,  int * nz) throw()
{
  // if (buffer) (*buffer) = (void * ) 
  // 		particle_data_->values(particle_descr_,particle_index_);
  // if (name)   (*name)   = particle_descr_->particle_name(particle_index_);
  // if (type) {

  //   precision_type precision = particle_descr_->precision(particle_index_);
  //   if (precision == precision_default) precision = default_precision;

  //   switch (precision) {
  //   case precision_single: (*type) = scalar_type_float; break;
  //   case precision_double: (*type) = scalar_type_double; break;
  //   default:
  //     ERROR2 ("IoParticleData",
  // 	      "Unsupported precision type %d for particle %s",
  // 	      precision, particle_descr_->particle_name(particle_index_).c_str());
  //   }
  // }
  // int nbx,nby,nbz;
  // particle_data_->size(&nbx,&nby,&nbz);

  // int ngx=0,ngy=0,ngz=0;

  // particle_descr_->ghost_depth(particle_index_,&ngx,&ngy,&ngz);

  // if (particle_data_->ghosts_allocated()) {

  //   if (nxd) (*nxd) = nbx + 2*ngx;
  //   if (nyd) (*nyd) = nby + 2*ngy;
  //   if (nzd) (*nzd) = nbz + 2*ngz;

  //   // Exclude ghosts when writing

  //   //    if (nx) (*nx) = nbx;
  //   //    if (ny) (*ny) = nby;
  //   //    if (nz) (*nz) = nbz;

  //   // Include ghosts when writing

  //    if (nx) (*nx) = nbx + 2*ngx;
  //    if (ny) (*ny) = nby + 2*ngy;
  //    if (nz) (*nz) = nbz + 2*ngz;

  // } else {

  //   if (nxd) (*nxd) = nbx;
  //   if (nyd) (*nyd) = nby;
  //   if (nzd) (*nzd) = nbz;

  //   if (nx) (*nx) = nbx;
  //   if (ny) (*ny) = nby;
  //   if (nz) (*nz) = nbz;

  // }
}
//----------------------------------------------------------------------
