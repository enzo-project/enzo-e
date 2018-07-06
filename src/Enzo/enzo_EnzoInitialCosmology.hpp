// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialCosmology.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-07-24
/// @brief    [\ref Enzo] Declaration of the EnzoInitialCosmology class

#ifndef ENZO_ENZO_INITIAL_COSMOLOGY_HPP
#define ENZO_ENZO_INITIAL_COSMOLOGY_HPP

class EnzoInitialCosmology : public Initial {

  /// @class    EnzoInitialCosmology
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoInitialCosmology(int cycle, double time,
		       double gamma,
		       double temperature) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialCosmology);

  /// CHARM++ migration constructor
  EnzoInitialCosmology(CkMigrateMessage *m) 
    : Initial (m),
      gamma_(0.0),
      temperature_(0.0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: update whenever attributes change

    TRACEPUP;

    Initial::pup(p);

    p | gamma_;
    p | temperature_;
  }
  
public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  (
   Block * block,
   const FieldDescr * field_descr,
   const ParticleDescr * particle_descr,
   const Hierarchy * hierarchy
   ) throw();

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Ideal gas law constant
  double gamma_;

  /// Initial temperature
  double temperature_;

};

#endif /* ENZO_ENZO_INITIAL_COSMOLOGY_HPP */

