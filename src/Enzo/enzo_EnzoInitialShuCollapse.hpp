// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.hpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     7 June 2022
/// @brief    Initializer for the Shu Collapse problem as described
///           in Federrath et al 2010, ApJ, 713, 269.

#ifndef ENZO_ENZO_INITIAL_SHU_COLLAPSE_HPP
#define ENZO_ENZO_INITIAL_SHU_COLLAPSE_HPP

class EnzoInitialShuCollapse : public Initial {

  /// @class    EnzoInitialShuCollapse
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the Shu Collapse problem (Shu 1977)

public: // interface

  /// Constructor
  EnzoInitialShuCollapse
  (int cycle, double time,
   const double center[3],
   const double drift_velocity[3],
   double truncation_radius,
   double nominal_sound_speed,
   double instability_parameter,
   double external_density,
   bool central_sink_exists,
   double central_sink_mass) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialShuCollapse);

  /// CHARM++ migration constructor
  EnzoInitialShuCollapse(CkMigrateMessage *m)
    : Initial (m),
      truncation_radius_(0.0),
      nominal_sound_speed_(0.0),
      instability_parameter_(0.0),
      central_sink_exists_(false),
      central_sink_mass_(0.0)
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

  /// Truncation radius - must be less than half the total domain width
  double truncation_radius_;

  /// Nominal uniform sound speed of the gas used to initialise the total specific energy
  /// In practice the actual sound speed will be different
  double nominal_sound_speed_;

  /// Instability parameter - sphere is gravitationally unstable if this is greater than 2.0
  /// Determines density profile.
  double instability_parameter_;

  /// Density outside of the truncation radius
  double external_density_;

  /// Is there are sink particle initialised at the centre of collapse?
  bool central_sink_exists_;

  /// Mass of the central sink particle. Ignored is central_sink_exists_ is false
  double central_sink_mass_;
};

#endif /* ENZO_ENZO_INITIAL_SHU_COLLAPSE_HPP */
