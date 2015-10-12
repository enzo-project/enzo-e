// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 11:48:15 PDT 2014
/// @brief    Implementation of the ParticleDescr class
///
/// ParticleDescr is used to describe the content of ParticleBlock's data
/// Separate classes are used to avoid redundant storage between ParticleBlocks,
/// which are designed to be extremely memory-efficient.
///
/// Particles have different types (e.g. tracer, dark matter, etc.),
/// and different types have different attributes (e.g. position,
/// velocity, mass, etc.)  ParticleDescr objects store which particle
/// types each ParticleBlock contains, and provides operations to
/// assist in "decoding" the data stored in ParticleBlocks.

#include "data.hpp"

//----------------------------------------------------------------------

ParticleDescr::ParticleDescr
(
 int rank) throw ()
  : rank_(rank),
    type_(),
    attribute_size_(),
    attribute_(),
    attribute_offset_(),
    groups_()
{
  // Initialize default attribute sizes
  attribute_size_.resize (num_particle_attribute);
  attribute_size_[particle_attribute_unknown]  = 0;
  attribute_size_[particle_attribute_position] = 4; // 10 bits per 3D axis
  attribute_size_[particle_attribute_velocity] = 8; // 20 bits per 3D axis
  attribute_size_[particle_attribute_mass]     = 4;
  attribute_size_[particle_attribute_id]       = 4;

}

//----------------------------------------------------------------------

void ParticleDescr::pup (PUP::er &p)
{
  TRACEPUP;

  p | rank_;
  p | type_;
  p | attribute_;  
  p | attribute_size_;
  p | attribute_offset_;
  p | groups_;
}

//======================================================================

void ParticleDescr::add_attribute (int index_type, int attribute) throw (){
  if (index_type >= (int)attribute_.size()) {
    attribute_.resize(index_type + 1);
  }
  if (index_type >= (int)attribute_offset_.size()) {
    attribute_offset_.resize(index_type + 1);
  }

  // number of attributes for the given type before adding this one
  int na = attribute_[index_type].size();

  // store the attribute type
  attribute_[index_type].push_back(attribute);

  // compute the index offset for this attribute

  int index;

  if (na == 0) {
    
    index = 1;

  } else {
    int sum = attribute_offset_[index_type][na-1];
    index = sum + attribute_size(attribute_[index_type][na-1]);
  }
  attribute_offset_[index_type].push_back(index);
}

//----------------------------------------------------------------------

void ParticleDescr::print () const
{
  printf ("rank = %d\n",rank_);
  for (size_t i=0; i<attribute_size_.size(); i++) {
    printf ("attribute_size_[%lu] = %d\n",i,attribute_size_[i]);
  }
  for (size_t i=0; i<type_.size(); i++) {
    printf ("type_[%lu] = %s\n",i,type_[i].c_str());
  }
  for (size_t i=0; i<attribute_.size(); i++) {
    for (size_t j=0; j<attribute_[i].size(); j++) {
      printf ("attribute_      [%lu][%lu] = %d\n",i,j,attribute_[i][j]);
    }
  }
  for (size_t i=0; i<attribute_offset_.size(); i++) {
    for (size_t j=0; j<attribute_offset_[i].size(); j++) {
      printf ("attribute_offset_[%lu][%lu] = %d\n",i,j,attribute_offset_[i][j]);
    }
  }
}
