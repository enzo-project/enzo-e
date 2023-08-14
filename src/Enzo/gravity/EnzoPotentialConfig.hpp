// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoPotentialConfig.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     June 6, 2023
/// @brief    [\ref Enzo] Declaration of EnzoPotentialConfigGalaxy and the
///           EnzoPotentialConfigPointMass classes.
///
/// Analytic potentials used by EnzoMethodBackgroundAcceleration.
///
/// At the moment, these are just glorified structs that store a bunch of
/// parameter values, that are not actually attached to the functions that do
/// the main work. That functionality is defined in a separate class because we
/// want to support compatability with cosmological simulations.
///
///   - If we didn't care about supporting cosmological simulations, we could
///     just convert the values of the parameters once (right after parsing
///     the parameters) and be done with it.
///   - Unfortunately, in cosmological simulations the code-units are
///     time-dependent
///   - Thus, we define the EnzoPotentialConfig class to store the parameters
///     in their default units, define a function to convert the parameters to
///     code units (and store the result in a new EnzoPotentialConfig instance)
///     and use that that new EnzoPotentialConfig instance to configure
///     the potential calculation
///
/// The way that the EnzoMethodBackgroundAcceleration is currently implemented
/// requires the definition of the functions that perform the actual potential
/// calculation (from a EnzoPotentialConfig instance) to be visible in the same
/// source file as the rest of the definition of
/// EnzoMethodBackgroundAcceleration.
///   - This could be accomplished by defining the functions inline, but at the
///     moment, that seems a little unnecssary.

#ifndef ENZO_GRAVITY_ENZO_POTENTIAL_CONFIG_HPP
#define ENZO_GRAVITY_ENZO_POTENTIAL_CONFIG_HPP

struct EnzoPotentialConfigGalaxy {

  /// @class    EnzoPotentialConfigGalaxy
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] tracks the parameters needed for a background
  ///           potential of a Galaxy model

  double DM_mass;
  double DM_mass_radius;
  double stellar_r;
  double stellar_z;
  double stellar_mass;
  double bulge_mass;
  double bulgeradius;
  std::array<double,3> amom;
  double rcore;

  static EnzoPotentialConfigGalaxy from_config(const EnzoConfig* enzo_config) {
    double DM_mass = enzo_config->method_background_acceleration_DM_mass;
    double DM_mass_radius = enzo_config->method_background_acceleration_DM_mass_radius;
    double stellar_r = enzo_config->method_background_acceleration_stellar_scale_height_r;
    double stellar_z = enzo_config->method_background_acceleration_stellar_scale_height_z;
    double stellar_mass = enzo_config->method_background_acceleration_stellar_mass;
    double bulge_mass = enzo_config->method_background_acceleration_bulge_mass;
    double bulgeradius = enzo_config->method_background_acceleration_bulge_radius;
    std::array<double,3> amom
      = {enzo_config->method_background_acceleration_angular_momentum[0],
         enzo_config->method_background_acceleration_angular_momentum[1],
         enzo_config->method_background_acceleration_angular_momentum[2]};
    double rcore = enzo_config->method_background_acceleration_core_radius;

    ASSERT1("GalaxyModel::GalaxyModel",
            "DM halo mass (=%e code_units) must be positive and specified in "
            "units of solar masses",
            DM_mass, (DM_mass > 0));

    return {DM_mass, DM_mass_radius, stellar_r, stellar_z,
            stellar_mass, bulge_mass, bulgeradius, amom, rcore};
  }

  void pup(PUP::er &p) {
    p | DM_mass;
    p | DM_mass_radius;
    p | stellar_r;
    p | stellar_z;
    p | stellar_mass;
    p | bulge_mass;
    p | bulgeradius;
    p | amom;
    p | rcore;
  }
};

//---------------------------------------------------------------------

struct EnzoPotentialConfigPointMass {
  double mass;
  double rcore;

  static EnzoPotentialConfigPointMass from_config(const EnzoConfig* enzo_config){
    double mass = enzo_config->method_background_acceleration_mass;
    double rcore = enzo_config->method_background_acceleration_core_radius;
    return {mass, rcore};
  }

  void pup(PUP::er &p) {
    p | mass;
    p | rcore;
  }
};

#endif /* ENZO_GRAVITY_ENZO_POTENTIAL_CONFIG_HPP */
