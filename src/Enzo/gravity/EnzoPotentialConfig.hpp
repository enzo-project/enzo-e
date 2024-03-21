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
  ///
  /// If we ever start using the idea of a Galaxy model across multiple Methods,
  /// it might make sense to somehow convert this into a Physics class


  /// dark matter mass (default units of solar masses)
  double DM_mass;
  /// dark matter mass radius (default units of kpc)
  double DM_mass_radius;
  /// disk stellar scale height r (default units of kpc)
  double stellar_r;
  /// disk stellar scale height z (default units of kpc)
  double stellar_z;
  /// disk stellar mass (default units of solar masses)
  double stellar_mass;
  /// bulge mass (default units of solar masses)
  double bulge_mass;
  /// bulgeradius (default units of kpc)
  double bulgeradius;
  /// angular momentum vector
  std::array<double,3> amom;
  /// core radius (default units of kpc)
  double rcore;

  static EnzoPotentialConfigGalaxy from_parameters(ParameterGroup p) {
    double DM_mass = p.value_float("DM_mass",-1.0);
    double DM_mass_radius = p.value_float("DM_mass_radius", 0.0);
    double stellar_r = p.value_float("stellar_scale_height_r", 1.0E-10);
    double stellar_z = p.value_float("stellar_scale_height_z", 1.0E-10);
    double stellar_mass = p.value_float("stellar_mass", 0.0);
    double bulge_mass = p.value_float("bulge_mass", 0.0);
    double bulgeradius = p.value_float("bulge_radius", 1.0E-10);
    std::array<double,3> amom;
    for (int axis = 0; axis < 3; axis++) {
      amom[axis] = p.list_value_float(axis,"angular_momentum",0);
    }
    double rcore = p.value_float("core_radius", 1.0E-10);

    ASSERT1("GalaxyModel::GalaxyModel",
            "DM halo mass (=%e code_units) must be positive and specified in "
            "units of solar masses",
            DM_mass, (DM_mass > 0));

    return {DM_mass, DM_mass_radius, stellar_r, stellar_z,
            stellar_mass, bulge_mass, bulgeradius, amom, rcore};
  }

  /// static method to convert from default units to code units
  static EnzoPotentialConfigGalaxy to_codeU
  (const EnzoPotentialConfigGalaxy& pack_dfltU, const Units* units) noexcept
  {
    return
      { pack_dfltU.DM_mass * enzo_constants::mass_solar / units->mass(),
        pack_dfltU.DM_mass_radius * enzo_constants::kpc_cm / units->length(),
        pack_dfltU.stellar_r * enzo_constants::kpc_cm / units->length(),
        pack_dfltU.stellar_z * enzo_constants::kpc_cm / units->length(),
        pack_dfltU.stellar_mass * enzo_constants::mass_solar / units->mass(),
        pack_dfltU.bulge_mass * enzo_constants::mass_solar / units->mass(),
        pack_dfltU.bulgeradius * enzo_constants::kpc_cm / units->length(),
        pack_dfltU.amom,
        pack_dfltU.rcore * enzo_constants::kpc_cm / units->length() };
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
  /// @class    EnzoPotentialConfigPointMass
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] tracks the parameters needed for a background
  ///           potential of an isolated point-mass


  /// point-source mass (default units of solar masses)
  double mass;
  /// core-radius (default units of cm)
  double rcore;

  static EnzoPotentialConfigPointMass from_parameters(ParameterGroup p) {
    return {p.value_float("mass",0.0),
            p.value_float("core_radius", 1.0E-10)};
  }

  /// static method to convert from default units to code units
  static EnzoPotentialConfigPointMass to_codeU
  (const EnzoPotentialConfigPointMass& pack_dfltU, const Units* units) noexcept
  {
    return
      { pack_dfltU.mass * enzo_constants::mass_solar / units->mass(),
        pack_dfltU.rcore / units->length() };
  }

  void pup(PUP::er &p) {
    p | mass;
    p | rcore;
  }
};

#endif /* ENZO_GRAVITY_ENZO_POTENTIAL_CONFIG_HPP */
