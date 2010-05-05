// $Id: data.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     data.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    Declaration of the Data class

#ifndef DATA_HPP
#define DATA_HPP

#include <vector>

#include "particles.hpp"
#include "field.hpp"

class Data {

  /// @class    Data
  /// @ingroup  Data
  /// @brief    Container class for Particles and Fields

public: // interface

  /// Initialize the Data object
  Data()
    : field_(),
      particles_()
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
  // Particles
  //----------------------------------------------------------------------

  /// Add particles
  void add_particles (Particles * particles)
  { particles_.push_back (particles); };

  /// Return the number of particless
  int particles_count ()
  { return particles_.size(); };

  /// Return the ith particles
  Particles * particles (int i)
  { 
    return (0 <= i && i < particles_count()) ? particles_.at(i) : 0;
  }

  /// Return the named particles
  Particles * particles (std::string name)
  {
    for (int i=0; i<particles_count(); i++) {
      if (particles(i)->name() == name) {
	return particles_.at(i); 
      };
    }
    return 0;
  }


private: // attributes

  std::vector<Field *>     field_;
  std::vector<Particles *> particles_;

};

#endif /* DATA_HPP */

