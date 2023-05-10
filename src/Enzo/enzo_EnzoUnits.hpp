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

    p | cosmology_; // PUP::able
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

  /// Return density units scaling factor (virtual)
  virtual double density() const
  {
    return (cosmology_ == NULL) ?
      Units::density() : cosmology_->density_units();
  }

  /// Return velocity units scaling factor (virtual)
  virtual double velocity() const
  {
    return (cosmology_ == NULL) ?
      Units::velocity() : cosmology_->velocity_units();
  }

  /// Returns the scaling factor to be divided by temperature to convert from
  /// units of Kelvin to units of specific internal energy
  ///
  /// The original Enzo refered to this quantity as "temperature units", but in
  /// Enzo-E we primarily track temperature in units of Kelvin
  double kelvin_per_energy_units() const
  {
    // The Units base-class can't compute the temperature code units for the
    // non-cosmology case because it involves physical constants that are not
    // defined in the Cello layer
    double vel_units = velocity();
    return (enzo_constants::mass_hydrogen * (vel_units * vel_units) /
	    enzo_constants::kboltz);
  }

  /// Return photon number density units scaling factor (derived)
  double photon_number_density() const
  {
    // The ratio of ionizing photons emitted per hydrogen atom at reionization is
    // just 1.5-3 (Bolton and Haehnelt 2007), so it makes sense to define our
    // cosmological photon number density units as proportional to
    // the cosmological baryon density units to keep the numbers of order unity in code units.
    return (cosmology_ == NULL) ?
      1.0 / Units::volume() : 0.76 * cosmology_->density_units() /
                              enzo_constants::mass_hydrogen;
  }
  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// EnzoCosmology units parameters
  EnzoPhysicsCosmology * cosmology_;

};

#endif /* ENZO_ENZO_UNITS_HPP */

