// $Id: particles.hpp 1261 2010-03-03 00:14:11Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARTICLES_PARTICLES_HPP
#define PARTICLES_PARTICLES_HPP

/// @file     particles_particles.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Declaration for the Particles class

class Particles {

  /// @class    Particles
  /// @ingroup  Particles
  /// @brief    Interface for the Particles class

public: // public

  /// Initialize a Particles object
  Particles()
    : name_()
  { };

  /// Set Particles name
  void set_name (std::string name)
  { name_ = name; };

  /// Get Particles name
  std::string name ()
  { return name_; };

private: // attributes

  /// String defining the field's name
  std::string name_;

};

#endif /* PARTICLES_PARTICLES_HPP */

