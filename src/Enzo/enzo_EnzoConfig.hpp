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
  EnzoConfig(CkMigrateMessage *m)
    : Config (m),
      ppm_density_floor(0.0),
      ppm_diffusion(false),
      ppm_dual_energy(false),
      ppm_dual_energy_eta_1(0.0),
      ppm_dual_energy_eta_2(0.0),
      ppm_flattening(0),
      ppm_minimum_pressure_support_parameter(0),
      ppm_number_density_floor(0.0),
      ppm_pressure_floor(0.0),
      ppm_pressure_free(false),
      ppm_steepening(false),
      ppm_temperature_floor(0.0),
      ppm_use_minimum_pressure_support(false),
      ppm_mol_weight(0.0),
      field_gamma(0.0),
      // Cosmology (NOT ACCESSED)
      physics_cosmology(false),
      physics_cosmology_comoving_box_size(0.0),
      physics_cosmology_hubble_constant_now(0.0),
      physics_cosmology_initial_redshift(0.0),
      physics_cosmology_max_expansion_rate(0.0),
      physics_cosmology_omega_lamda_now(0.0),
      physics_cosmology_omega_matter_now(0.0),
      // EnzoInitialSedovArray[23]
      initial_sedov_rank(0),
      initial_sedov_radius_relative(0.0),
      initial_sedov_pressure_in(0.0),
      initial_sedov_pressure_out(0.0),
      initial_sedov_density(0.0),
      initial_turbulence_density(0.0),
      initial_turbulence_pressure(0.0),
      initial_turbulence_temperature(0.0),
      // EnzoInitialPm
      initial_pm_field(""),
      initial_pm_mpp(0.0),
      // EnzoInitialSoup
      initial_soup_rank(0),
      initial_soup_file(""),
      initial_soup_pressure_in(0.0),
      initial_soup_pressure_out(0.0),
      initial_soup_density(0.0),
      // EnzoProlong
      interpolation_method(""),
      // EnzoMethodHeat
      method_heat_alpha(0.0),
      // EnzoMethodNull
      method_null_dt(0.0),
      // EnzoMethodTurbulence
      method_turbulence_edot(0.0),
      method_turbulence_mach_number(0.0),
      // EnzoMethodGravity
      method_gravity_grav_const(0.0),
      method_gravity_solver(),
      // EnzoMethodGravityCg
      method_gravity_cg_iter_max(0),
      method_gravity_cg_res_tol(0.0),
      method_gravity_cg_grav_const(0.0),
      method_gravity_cg_diag_precon(false),
      method_gravity_cg_monitor_iter(0),
      // EnzoMethodGravityBiCGStab
      method_gravity_bicgstab_iter_max(0),
      method_gravity_bicgstab_res_tol(0.0),
      method_gravity_bicgstab_grav_const(0.0),
      method_gravity_bicgstab_diag_precon(false),
      method_gravity_bicgstab_monitor_iter(0),
      // EnzoMethodGravityMlat
      // EnzoMethodGravityMg0
      method_gravity_mg_type(""),
      method_gravity_mg_iter_max(0),
      method_gravity_mg_res_tol(0.0),
      method_gravity_mg_grav_const(0.0),
      method_gravity_mg_monitor_iter(0),
      method_gravity_mg_smooth(""),
      method_gravity_mg_smooth_weight(0.0),
      method_gravity_mg_smooth_pre(0),
      method_gravity_mg_smooth_coarse(0),
      method_gravity_mg_smooth_post(0),
      method_gravity_mg_restrict(""),
      method_gravity_mg_prolong(""),
      method_gravity_mg_min_level(0),
      method_gravity_mg_max_level(0),
      // EnzoMethodPm
      method_pm_deposit_type(""),
      method_pm_update_max_dt(0.0)
  {
    for (int axis=0; axis<3; axis++) {
      initial_sedov_array[axis] = 0;
      initial_soup_array[axis] = 0;
      initial_soup_d_pos[axis] = 0;
      initial_soup_d_size[axis] = 0;
    }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();

  /// Write configuration to a file
  void write (FILE * fp);
  
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
  int                        initial_sedov_rank;
  int                        initial_sedov_array[3];
  double                     initial_sedov_radius_relative;
  double                     initial_sedov_pressure_in;
  double                     initial_sedov_pressure_out;
  double                     initial_sedov_density;

  double                     initial_turbulence_density;
  double                     initial_turbulence_pressure;
  double                     initial_turbulence_temperature;

  std::string                initial_pm_field;
  double                     initial_pm_mpp;
  int                        initial_pm_level;

  int                        initial_soup_rank;
  std::string                initial_soup_file;
  bool                       initial_soup_rotate;
  int                        initial_soup_array[3];
  double                     initial_soup_d_pos[3];
  double                     initial_soup_d_size[3];
  double                     initial_soup_pressure_in;
  double                     initial_soup_pressure_out;
  double                     initial_soup_density;

  // EnzoProlong
  std::string                interpolation_method;

  // EnzoMethodHeat
  double                     method_heat_alpha;

  // EnzoMethodNull
  double                     method_null_dt;

  // EnzoMethodTurbulence
  double                     method_turbulence_edot;
  double                     method_turbulence_mach_number;

  // EnzoMethodGravity
  double                     method_gravity_grav_const;
  std::string                method_gravity_solver;

  // EnzoMethodGravityCg
  int                        method_gravity_cg_iter_max;
  double                     method_gravity_cg_res_tol;
  double                     method_gravity_cg_grav_const;
  bool                       method_gravity_cg_diag_precon;
  int                        method_gravity_cg_monitor_iter;

  // EnzoMethodGravityBiCGStab
  int                        method_gravity_bicgstab_iter_max;
  double                     method_gravity_bicgstab_res_tol;
  double                     method_gravity_bicgstab_grav_const;
  bool                       method_gravity_bicgstab_diag_precon;
  int                        method_gravity_bicgstab_monitor_iter;

  // EnzoMethodGravityMlat
  // EnzoMethodGravityMg0
  std::string                method_gravity_mg_type;
  int                        method_gravity_mg_iter_max;
  double                     method_gravity_mg_res_tol;
  double                     method_gravity_mg_grav_const;
  int                        method_gravity_mg_monitor_iter;
  std::string                method_gravity_mg_smooth;
  double                     method_gravity_mg_smooth_weight;
  int                        method_gravity_mg_smooth_pre;
  int                        method_gravity_mg_smooth_coarse;
  int                        method_gravity_mg_smooth_post;
  std::string                method_gravity_mg_restrict;
  std::string                method_gravity_mg_prolong;
  int                        method_gravity_mg_min_level;
  int                        method_gravity_mg_max_level;

  // EnzoMethodPm

  std::string                method_pm_deposit_type;
  double                     method_pm_update_max_dt;

#ifdef CONFIG_USE_GRACKLE

  // EnzoMethodGrackle

  code_units      method_grackle_units;
  chemistry_data  method_grackle_chemistry;

#endif /* CONFIG_USE_GRACKLE */

};

extern EnzoConfig g_enzo_config;

#endif /* PARAMETERS_ENZO_CONFIG_HPP */

