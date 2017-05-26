// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoUnits.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-05-22
/// @brief    [\ref Enzo] Declaration of the EnzoUnits class

#ifndef ENZO_ENZO_UNITS_HPP
#define ENZO_ENZO_UNITS_HPP

class EnzoUnits {

  /// @class    EnzoUnits
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoUnits() throw()
  : time_(1.0),
    mass_(1.0),
    length_(1.0),
    density_(1.0),
    temperature_(1.0),
    velocity_(1.0),
    cosmology_(NULL)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) {
    TRACEPUP;

    p | time_;
    p | mass_;
    p | length_;
    p | density_;
    p | temperature_;
    p | velocity_;

    bool not_null = (cosmology_ != NULL);
    p | not_null;
    if (not_null) {
      if (p.isUnpacking()) cosmology_ = new EnzoCosmology;
      p | *cosmology_;
    } else {
      cosmology_ = NULL;
    }

  }

  double time() const
  { return time_; }
  double mass() const
  { return mass_; }
  double length() const
  { return length_; }
  double density() const
  { return density_; }
  double temperature() const
  { return temperature_; }
  double velocity() const
  { return length_ / time_; }
  
  /// Set units given mass
  void set_given_mass (double length, double mass, double time)
  {
    length_ = length;
    mass_   = mass;
    time_   = time;

    density_ = mass / length / length / length;
  }

  /// Set units given density
  void set_given_density (double length, double density, double time)
  {
    length_  = length;
    density_ = density;
    time_    = time;

    mass_ = density * length * length * length;
  }

  void set_cgs();

  /// Update current units for cosmological problems
  void update_cosmology (double time);
  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  double time_;
  double mass_;
  double length_;
  double density_;
  double temperature_;
  double velocity_;

  /// EnzoCosmology units parameters
  EnzoCosmology * cosmology_;
  
};

#endif /* ENZO_ENZO_UNITS_HPP */

