// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoConfig.hpp
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
  double                     ppm_temperature_floor;
  bool                       ppm_use_minimum_pressure_support;
  double                     ppm_mol_weight;

  double                     field_gamma;

  // Cosmology (NOT ACCESSED)
  bool                       physics_cosmology;
  double                     physics_cosmology_comoving_box_size;
  double                     physics_cosmology_hubble_constant_now;
  double                     physics_cosmology_initial_redshift;
  double                     physics_cosmology_max_expansion_rate;
  double                     physics_cosmology_omega_lamda_now;
  double                     physics_cosmology_omega_matter_now;

  // EnzoInitialSedovArray[23]
  int                        initial_sedov_array[3];
  double                     initial_sedov_radius_relative;
  double                     initial_sedov_pressure_in;
  double                     initial_sedov_pressure_out;
  double                     initial_sedov_density;

  double                     initial_turbulence_density;
  double                     initial_turbulence_pressure;
  double                     initial_turbulence_temperature;

  // EnzoProlong
  std::string                interpolation_method;

  // EnzoMethodHeat
  double                     method_heat_alpha;

  // EnzoMethodNull
  double                     method_null_dt;

  // EnzoMethodTurbulence
  double                     method_turbulence_edot;
  double                     method_turbulence_mach_number;

  // EnzoMethodGravityCg
  int                        method_gravity_cg_iter_max;
  double                     method_gravity_cg_res_tol;

#ifdef CONFIG_USE_GRACKLE

  // EnzoMethodGrackle

  code_units      method_grackle_units;
  chemistry_data  method_grackle_chemistry;

#endif /* CONFIG_USE_GRACKLE */

};

#endif /* PARAMETERS_ENZO_CONFIG_HPP */

