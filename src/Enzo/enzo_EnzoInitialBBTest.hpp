// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialBBTest.hpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     10 June 2022
/// @brief    Initializer for the "BB Test" problem as described
///           in Federrath et al 2010, ApJ, 713, 269.

#ifndef ENZO_ENZO_INITIAL_BB_TEST_HPP
#define ENZO_ENZO_INITIAL_BB_TEST_HPP

class EnzoInitialBBTest : public Initial {

  /// @class    EnzoInitialBBTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the Shu Collapse problem (Shu 1977)

public: // interface

  /// Constructor
  EnzoInitialBBTest
  (int cycle, double time,
   const double center[3],
   const double drift_velocity[3],
   double mean_density,
   double fluctuation_amplitude,
   double truncation_radius,
   double nominal_sound_speed,
   double angular_rotation_velocity,
   double external_density) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialBBTest);

  /// CHARM++ migration constructor
  EnzoInitialBBTest(CkMigrateMessage *m)
    : Initial (m),
      truncation_radius_(0.0),
      nominal_sound_speed_(0.0),
      mean_density_(0.0),
      fluctuation_amplitude_(0.0),
      angular_rotation_velocity_(0.0),
      external_density_(0.0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

  private: // attributes

  /// Location of the centre of collapse
  double center_[3];

  /// Drift velocity of the whole system
  double drift_velocity_[3];

  /// Mean initial density
  double mean_density_;

  /// The amplitude of the density fluctuation
  double fluctuation_amplitude_;

  /// Truncation radius - must be less than half the total domain width
  double truncation_radius_;

  /// Nominal uniform sound speed of the gas used to initialise the total specific energy
  /// In practice the actual sound speed will be different
  double nominal_sound_speed_;

  /// The angular rotation velocity of the sphere (units are radians per unit time)
  double angular_rotation_velocity_;

  /// The density outside of the truncation radius;
  double external_density_;
};

#endif /* ENZO_ENZO_INITIAL_BB_TEST_HPP */
