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

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoConfig);

  /// CHARM++ migration constructor
  EnzoConfig(CkMigrateMessage *m) : Config (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();

public: // attributes

  // NOTE: change pup() function whenever attributes change

  // EnzoMethodPpm

  double                     ppm_density_floor;
  bool                       ppm_diffusion;
  bool                       ppm_dual_energy;
  double                     ppm_dual_energy_eta_1;
  double                     ppm_dual_energy_eta_2;
  int                        ppm_flattening;
  int                        ppm_minimum_pressure_support_parameter;
  double                     ppm_number_density_floor;
  double                     ppm_pressure_floor;
  bool                       ppm_pressure_free;
  bool                       ppm_steepening;
  float                      ppm_temperature_floor;
  bool                       ppm_use_minimum_pressure_support;

  double                     field_gamma;

  // Cosmology (NOT ACCESSED)
  bool                       cosmology;
  double                     cosmology_comoving_box_size;
  double                     cosmology_hubble_constant_now;
  double                     cosmology_initial_redshift;
  double                     cosmology_max_expansion_rate;
  double                     cosmology_omega_lamda_now;
  double                     cosmology_omega_matter_now;

  // EnzoInitialSedovArray[23]
  int                        sedov_array[3];
  double                     sedov_radius_relative;
  double                     sedov_pressure_in;
  double                     sedov_pressure_out;
  double                     sedov_density;

  // EnzoProlong
  std::string                interpolation_method;

  // EnzoMethodHeat
  double                     method_heat_alpha;
};

#endif /* PARAMETERS_ENZO_CONFIG_HPP */

