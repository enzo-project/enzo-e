// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoUnits.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-05-22
/// @brief    [\ref Enzo] Declaration of the EnzoUnits class

#ifndef ENZO_ENZO_UNITS_HPP
#define ENZO_ENZO_UNITS_HPP

class EnzoUnits : public Units {

  /// @class    EnzoUnits
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoUnits() throw()
  : Units(),
    cosmology_(NULL)
  {
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoUnits);
  
  /// CHARM++ migration constructor for PUP::able
  EnzoUnits (CkMigrateMessage *m)
    : Units (m),
      cosmology_(NULL)
  {  }

  /// Virtual destructor
  virtual ~EnzoUnits ()
  { }
  
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) {

    TRACEPUP;

    Units::pup(p);

    bool not_null = (cosmology_ != NULL);
    p | not_null;
    if (not_null) {
      if (p.isUnpacking()) cosmology_ = new EnzoPhysicsCosmology;
      p | *cosmology_;
    } else {
      cosmology_ = NULL;
    }

  }

  /// Update current units for cosmological problems
  void set_current_time (double time)
  {
    if (cosmology_ != NULL) {
      cosmology_->set_current_time(time);
    }
  }

  /// Return the Cosmology object (may be NULL)
  EnzoPhysicsCosmology * cosmology()
  { return cosmology_; }
  
  /// Initialize the cosmology object, if any
  void set_cosmology(EnzoPhysicsCosmology * cosmology)
  { cosmology_ = cosmology; }
  
  /// Return magnetic units scaling factor
  inline double magnetic() const
  { return std::sqrt(pressure()*4.0*(cello::pi)); }

public: // virtual methods
 
  /// Return time units scaling factor (virtual)
  virtual double time() const
  {
    return (cosmology_ == NULL) ?
      Units::time() : cosmology_->time_units();
  }
  
  /// Return mass units scaling factor (virtual)
  virtual double mass() const
  {
    return (cosmology_ == NULL) ?
      Units::mass() : cosmology_->mass_units();
  }

  /// Return length units scaling factor (virtual)
  virtual double length() const
  {
    return (cosmology_ == NULL) ?
      Units::length() : cosmology_->length_units();
  }

  /// Return temperature units scaling factor (virtual)
  virtual double temperature() const
  {
    return (cosmology_ == NULL) ?
      Units::temperature() : cosmology_->temperature_units();
  }
  
  /// Return velocity units scaling factor (virtual)
  virtual double velocity() const
  {
    return (cosmology_ == NULL) ?
      Units::velocity() : cosmology_->velocity_units();
  }
  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// EnzoCosmology units parameters
  EnzoPhysicsCosmology * cosmology_;

};

#endif /* ENZO_ENZO_UNITS_HPP */

