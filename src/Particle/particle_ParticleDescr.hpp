// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARTICLE_PARTICLE_DESCR_HPP
#define PARTICLE_PARTICLE_DESCR_HPP

/// @file     particle_ParticleDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Declaration for the ParticleDescr class

class ParticleDescr {

  /// @class    ParticleDescr
  /// @ingroup  Particle
  /// @brief    Interface for the ParticleDescr class

public: // public

  /// Initialize a ParticleDescr object
  ParticleDescr()
    : name_()
  { };

  /// Set ParticleDescr name
  void set_name (std::string name)
  { name_ = name; };

  /// Get ParticleDescr name
  std::string name ()
  { return name_; };

private: // attributes

  /// String defining the Particle's name
  std::string name_;

};

#endif /* PARTICLE_PARTICLE_DESCR_HPP */

