// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_EnzoConfig.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-02
/// @brief    [\ref Parameters] Declaration of the EnzoConfig class
///

#ifndef PARAMETERS_ENZO_CONFIG_HPP
#define PARAMETERS_ENZO_CONFIG_HPP

#define MAX_FIELDS      30
#define MAX_FILE_GROUPS 10

class Parameters;

class EnzoConfig : public Config {

  /// @class    EnzoConfig
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Declaration of Enzo configuration class

public: // interface

  /// Constructor
  EnzoConfig() throw();

  /// Destructor
  ~EnzoConfig() throw();

  /// Copy constructor
  EnzoConfig(const EnzoConfig & config) throw();

  /// Assignment operator
  EnzoConfig & operator= (const EnzoConfig & config) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
#endif

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();

public: // attributes

  // NOTE: change pup() function whenever attributes change

  double                     enzo_ppm_density_floor;
  bool                       enzo_ppm_diffusion;
  bool                       enzo_ppm_dual_energy;
  double                     enzo_ppm_dual_energy_eta_1;
  double                     enzo_ppm_dual_energy_eta_2;
  int                        enzo_ppm_flattening;
  int                        enzo_ppm_minimum_pressure_support_parameter;
  double                     enzo_ppm_number_density_floor;
  double                     enzo_ppm_pressure_floor;
  bool                       enzo_ppm_pressure_free;
  bool                       enzo_ppm_steepening;
  float                      enzo_ppm_temperature_floor;
  bool                       enzo_ppm_use_minimum_pressure_support;

  bool                       physics_cosmology;
  double                     physics_cosmology_comoving_box_size;
  double                     physics_cosmology_hubble_constant_now;
  double                     physics_cosmology_initial_redshift;
  double                     physics_cosmology_max_expansion_rate;
  double                     physics_cosmology_omega_lamda_now;
  double                     physics_cosmology_omega_matter_now;
  double                     physics_gamma;


};

#endif /* PARAMETERS_ENZO_CONFIG_HPP */

