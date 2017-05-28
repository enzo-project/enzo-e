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
    : Units (m)
  {  }

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
  void update_cosmology (double time);

  /// Initialize the cosmology object, if any
  void set_cosmology(Physics * cosmology)
  { cosmology_ = cosmology; }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// EnzoCosmology units parameters
  Physics * cosmology_;
  
};

#endif /* ENZO_ENZO_UNITS_HPP */

