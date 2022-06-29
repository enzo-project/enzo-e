// See LICENSE_CELLO file for license and copyright information

/// @file     enzo.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-06-14
/// @brief    Physical constants

#ifndef ENZO_CONSTANTS_HPP
#define ENZO_CONSTANTS_HPP

/// Namespace for global constants and functions
namespace enzo_constants {

  // ergs per eV
  const double erg_eV = 1.60217653E-12;

  // eV per erg
  const double eV_erg = 6.24150948E11;

  // Boltzman constant in CGS
  const double kboltz = 1.3806504e-16;

  // Solar mass in CGS
  const double mass_solar = 1.98841586e33;

  // Solar metallicity
  const double metallicity_solar = 0.012;

  // Hydrogen mass in CGS
  const double mass_hydrogen = 1.67262171e-24;

  // Electron mass in CGS
  const double mass_electron = 9.10938291E-28;

  // parsec in CGS
  const double pc_cm  = 3.0856775809623245E18;

  // Kiloparsec in CGS
  const double kpc_cm = 3.0856775809623245E21;

  // Megaparsec in CGS
  const double Mpc_cm = 3.0856775809623245E24;

  // speed of light in CGS
  const double clight = 29979245800.0;

  // Gravitational constant in CGS
  // Note: in non-cosmological simulations, the Gravity solver uses the value specified by 
  // Method:gravity:grav_constant
  const double grav_constant = 6.67384E-8;

  // year in seconds
  const double yr_s = 3.1556952E7;

  // kyr in seconds
  const double kyr_s = 3.1556952E10;

  // Myr in seconds
  const double Myr_s = 3.1556952E13;

  // Approximate mean molecular weight of metals
  const double mu_metal = 16.0;

  // 100 (km / s) / Mpc in inverse seconds
  const double H0_over_h = 1.0e7 / Mpc_cm;

};

#endif /* ENZO_CONSTANTS_HPP */
