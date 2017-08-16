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
    length_(1.0)
  {
  }

  /// Virtual destructor
  virtual ~Units ()
  { }
  
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

  }

  /// Set units using mass
  void set_using_mass (double length, double mass, double time)
  {
    length_ = length;
    mass_   = mass;
    time_   = time;
  }

  /// Set units using density
  void set_using_density (double length, double density, double time)
  {
    length_  = length;
    mass_    = density * length * length * length;
    time_    = time;
  }

  /// Return volume units scaling factor (derived)
  inline double volume() const
  { return length()*length()*length(); }

  /// Return density units scaling factor (derived)
  inline double density() const
  { return mass() / volume(); }

  /// Return acceleration units scaling factor (derived)
  inline double acceleration() const
  { return length() / time() / time(); }

  /// Return pressure units scaling factor (derived)
  inline double pressure() const
  { return density() * velocity() * velocity(); }

public: // virtual methods

  /// Return time units scaling factor (virtual)
  virtual double time() const
  { return time_; }

  /// Return mass units scaling factor (virtual)
  virtual double mass() const
  { return mass_; }

  /// Return length units scaling factor (virtual)
  virtual double length() const
  { return length_; }
  
  /// Return temperature units scaling factor (derived)
  virtual double temperature() const
  { return (cello::mass_hydrogen)*std::pow(length()/time(),2)/(cello::k); }

  /// Return velocity units scaling factor (derived)
  virtual double velocity() const
  { return length() / time(); }
  
private: // attributes

  // NOTE: change pup() function whenever attributes change

  double time_;
  double mass_;
  double length_;

};

#endif /* PROBLEM_UNITS_HPP */
