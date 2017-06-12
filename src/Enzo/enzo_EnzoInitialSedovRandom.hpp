// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialSedovRandom.hpp
/// @author   Thomas Bolden (boldenth@gmail.com)
/// @date     Sat Mar  25 14:30:00 EST 2017
/// @brief    [\ref Enzo] 3D array of Randomized Sedov Blasts initial conditions

#ifndef ENZO_ENZO_INITIAL_SEDOV_RANDOM_HPP
#define ENZO_ENZO_INITIAL_SEDOV_RANDOM_HPP

class EnzoInitialSedovRandom : public Initial {

  /// @class    EnzoInitialSedovRandom
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] initial conditions for 3D array of Randomized Sedov Blasts 

public: // interface

  /// CHARM++ constructor
  EnzoInitialSedovRandom() throw()
  : Initial (),
    half_empty_(false),
    grackle_cooling_(false), // not implimented
    max_blasts_(0),
    radius_relative_(0.0),
    pressure_in_(0.0),
    pressure_out_(0.0),
    density_(0.0),
    te_multiplier_(0)
  { }
  
  /// Constructor
  EnzoInitialSedovRandom(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialSedovRandom);

  /// CHARM++ migration constructor
  EnzoInitialSedovRandom(CkMigrateMessage *m)
    : Initial (m),
      half_empty_(false),
      grackle_cooling_(false), // not implimented
      max_blasts_(0),
      radius_relative_(0.0),
      pressure_in_(0.0),
      pressure_out_(0.0),
      density_(0.0),
      te_multiplier_(0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( Block * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy * hierarchy ) throw();

private: // attributes

  /// Size of the array of Sedov blasts
  int array_[3];

  /// Whether to unbalance load with some boxes empty
  bool half_empty_;

  /// Whether to use Grackle for chemistry and radiative cooling
  /// NOT IMPLIMENTED
  bool grackle_cooling_;

  /// Maximum number of explosions per box
  int max_blasts_;

  /// Relative radius
  double radius_relative_;

  /// Internal and external pressure
  double pressure_in_;
  double pressure_out_;

  /// Initial density
  double density_;

  /// Variation of energy in sedov explosions
  int te_multiplier_;

};

#endif /* ENZO_ENZO_INITIAL_SEDOV_RANDOM_HPP */

