// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Units.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2017-05-26
/// @brief    [\ref Problem] Declaration for the Units component

#ifndef PROBLEM_UNITS_HPP
#define PROBLEM_UNITS_HPP

class Units : public PUP::able {

  /// @class    EnzoUnits
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  Units() throw()
  : time_(1.0),
    mass_(1.0),
    length_(1.0),
    density_(1.0),
    temperature_(1.0),
    velocity_(1.0)
  {
    CkPrintf ("TRACE Units\n");
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(Units);

  /// CHARM++ migration constructor for PUP::able
  Units (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) {

    PUP::able::pup(p);

    TRACEPUP;
    
    p | time_;
    p | mass_;
    p | length_;
    p | density_;
    p | temperature_;
    p | velocity_;

  }

  virtual double time() const        { return time_; }
  virtual double mass() const        { return mass_; }
  virtual double length() const      { return length_; }
  virtual double density() const     { return density_; }
  virtual double temperature() const { return temperature_; }
  virtual double velocity() const    { return length_ / time_; }
  
  /// Set units using mass
  void set_using_mass (double length, double mass, double time)
  {
    length_ = length;
    mass_   = mass;
    time_   = time;

    density_ = mass / length / length / length;
  }

  /// Set units using density
  void set_using_density (double length, double density, double time)
  {
    length_  = length;
    density_ = density;
    time_    = time;

    mass_ = density * length * length * length;
  }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  double time_;
  double mass_;
  double length_;
  double density_;
  double temperature_;
  double velocity_;

};

#endif /* PROBLEM_UNITS_HPP */
