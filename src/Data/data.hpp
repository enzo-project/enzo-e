// $Id: data.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef DATA_HPP
#define DATA_HPP

#include "particles.hpp"
#include "field.hpp"

/// @file     data.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    Declaration of the Data class

class Data {

  /// @class    Data
  /// @ingroup  Data
  /// @brief    Container class for Particles and Fields

public: // interface

  /// Initialize the Data object
  Data()
    : field_(0),
      particles_(0)
  { };

  /// Delete the Data object
  ~Data() { };

  /// Set fields
  void set_fields (Field * field)
  { field_ = field; };

  /// Get fields
  void get_fields (Field ** field)
  { *field = field_; };

  /// Set particles
  void set_particles (Particles * particles)
  { particles_ = particles; };

  /// Get particles
  void get_particles (Particles ** particles)
  { *particles = particles_; };

private: // attributes

  Field     * field_;
  Particles * particles_;

};

#endif /* DATA_HPP */

