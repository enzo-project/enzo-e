// See LICENSE_CELLO file for license and copyright information

/// @file     data_Particle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-03-31
/// @brief    

#include "data.hpp"
#include "charm_simulation.hpp"

//----------------------------------------------------------------------

int Particle::insert_particles (int it, int np)
{
  particle_data_->insert_particles (particle_descr_, it, np);
  
  return np;
}

//----------------------------------------------------------------------

int Particle::delete_particles (int it, int ib, const bool * mask)
{
  int np = particle_data_->delete_particles (particle_descr_,it,ib,mask); 

  return np;
}

//----------------------------------------------------------------------

Particle::~Particle() throw()
{
}
