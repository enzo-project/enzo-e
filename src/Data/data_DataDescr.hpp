// $Id: data_data.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     data_data.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    Declaration of the DataDescr class

#ifndef DATA_DATA_HPP
#define DATA_DATA_HPP

#include <vector>

#include "particle.hpp"
#include "field.hpp"

class DataDescr {

  /// @class    DataDescr
  /// @ingroup  Data
  /// @brief    Container class for Particles and Fields

public: // interface

  /// Initialize the DataDescr object
  DataDescr()
    : field_(),
      particle_()
  { };

  //----------------------------------------------------------------------
  // Fields
  //----------------------------------------------------------------------

  /// Add field
  void add_field (Field * field)
  { field_.push_back (field); };

  /// Return the number of fields
  int field_count ()
  { return field_.size(); };

  /// Return the ith field
  Field * field (int i)
  {
    return (0 <= i && i < field_count()) ? field_.at(i) : 0;
  }

  /// Return the named field
  Field * field (std::string name)
  {
    for (int i=0; i<field_count(); i++) {
      if (field(i)->name() == name) {
	return field_.at(i); 
      };
    }
    return 0;
  }

  //----------------------------------------------------------------------
  // Particle
  //----------------------------------------------------------------------

  /// Add particle
  void add_particle (ParticleDescr * particle)
  { particle_.push_back (particle); };

  /// Return the number of particles
  int particle_count ()
  { return particle_.size(); };

  /// Return the ith particle
  ParticleDescr * particle (int i)
  { 
    return (0 <= i && i < particle_count()) ? particle_.at(i) : 0;
  }

  /// Return the named particle
  ParticleDescr * particle (std::string name)
  {
    for (int i=0; i<particle_count(); i++) {
      if (particle(i)->name() == name) {
	return particle_.at(i); 
      };
    }
    return 0;
  }


private: // attributes

  std::vector<Field *>    field_;
  std::vector<ParticleDescr *> particle_;

};

#endif /* DATA_DATA_HPP */

