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
  virtual ~EnzoConfig() throw();

  /// Copy constructor
  EnzoConfig(const EnzoConfig & config) throw();

  /// Assignment operator
  EnzoConfig & operator= (const EnzoConfig & config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoConfig);

  /// CHARM++ migration constructor
  EnzoConfig(CkMigrateMessage *m)
    : Config (m),
      adapt_mass_type(),
      ppm_diffusion(0),
      ppm_dual_energy(false),
      ppm_dual_energy_eta_1(0.0),
      ppm_dual_energy_eta_2(0.0),
      ppm_flattening(0),
      ppm_minimum_pressure_support_parameter(0),
      ppm_number_density_floor(0.0),
      ppm_density_floor(0.0),
      ppm_pressure_floor(0.0),
      ppm_pressure_free(false),
      ppm_temperature_floor(0.0),
      ppm_steepening(false),
      ppm_use_minimum_pressure_support(false),
      ppm_mol_weight(0.0),
      field_gamma(0.0),
      // Cosmology
      physics_cosmology(false),
      physics_cosmology_hubble_constant_now(0.0),
      physics_cosmology_omega_matter_now(0.0),
      physics_cosmology_omega_lamda_now(0.0),
      physics_cosmology_omega_baryon_now(1.0),
      physics_cosmology_omega_cdm_now(0.0),
      physics_cosmology_comoving_box_size(0.0),
      physics_cosmology_max_expansion_rate(0.0),
      physics_cosmology_initial_redshift(0.0),
      physics_cosmology_final_redshift(0.0),
      physics_gravity(false),
      // EnzoInitialCosmology
      initial_cosmology_temperature(0.0),
      // EnzoInitialCollapse
      initial_collapse_rank(0),
      initial_collapse_radius_relative(0.0),
      initial_collapse_particle_ratio(0.0),
      initial_collapse_mass(0.0),
      initial_collapse_temperature(0.0),
      // EnzoInitialMusic
      initial_music_field_files(),
      initial_music_field_datasets(),
      initial_music_field_names(),
      initial_music_field_coords(),
      initial_music_particle_files(),
      initial_music_particle_datasets(),
      initial_music_particle_coords(),
      initial_music_particle_types(),
      initial_music_particle_attributes(),
      // EnzoInitialPm
      initial_pm_field(""),
      initial_pm_mpp(0.0),
      initial_pm_level(0),
      // EnzoInitialSedovArray[23]
      initial_sedov_rank(0),
      initial_sedov_radius_relative(0.0),
      initial_sedov_pressure_in(0.0),
      initial_sedov_pressure_out(0.0),
      initial_sedov_density(0.0),
      // EnzoInitialSedovRandom
      initial_sedov_random_half_empty(false),
      initial_sedov_random_grackle_cooling(false),
      initial_sedov_random_max_blasts(0),
      initial_sedov_random_radius_relative(0.0),
      initial_sedov_random_pressure_in(0.0),
      initial_sedov_random_pressure_out(0.0),
      initial_sedov_random_density(0.0),
      initial_sedov_random_te_multiplier(0),
      // EnzoInitialSoup
      initial_soup_rank(0),
      initial_soup_file(""),
      initial_soup_rotate(false),
      initial_soup_pressure_in(0.0),
      initial_soup_pressure_out(0.0),
      initial_soup_density(0.0),
      // EnzoInitialTurbulence
      initial_turbulence_density(0.0),
      initial_turbulence_pressure(0.0),
      initial_turbulence_temperature(0.0),
      // EnzoProlong
      interpolation_method(""),
      // EnzoMethodHeat
      method_heat_alpha(0.0),
      // EnzoMethodHydro
      method_hydro_method(""),
      method_hydro_dual_energy(false),
      method_hydro_dual_energy_eta_1(0.0),
      method_hydro_dual_energy_eta_2(0.0),
      method_hydro_reconstruct_method(""),
      method_hydro_reconstruct_conservative(false),
      method_hydro_reconstruct_positive(false),
      method_hydro_riemann_solver(""),
      // EnzoMethodNull
      method_null_dt(0.0),
      // EnzoMethodTurbulence
      method_turbulence_edot(0.0),
      method_turbulence_mach_number(0.0),
      // EnzoMethodGravity
      method_gravity_grav_const(0.0),
      method_gravity_solver(""),
      method_gravity_order(4),
      method_gravity_accumulate(false),
      // EnzoMethodPmDeposit
      method_pm_deposit_alpha(0.5),
      // EnzoMethodPmUpdate
      method_pm_update_max_dt(0.0),
      // EnzoSolverMg0
      solver_pre_smooth(),
      solver_post_smooth(),
      solver_last_smooth(),
      solver_coarse_solve(),
      solver_domain_solve(),
      solver_weight(),
      solver_restart_cycle(),
      // EnzoSolver<Krylov>
      solver_precondition(),
      solver_local(),
      solver_coarse_level(),
      solver_is_unigrid(),
      // EnzoStopping
      stopping_redshift()
      
  {
    for (int axis=0; axis<3; axis++) {
      initial_sedov_array[axis] = 0;
      initial_sedov_random_array[axis] = 0;
      initial_soup_array[axis] = 0;
      initial_soup_d_pos[axis] = 0;
      initial_soup_d_size[axis] = 0;
      initial_collapse_array[axis] = 0;
    }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();
  
public: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Refine

  std::vector <std::string>  adapt_mass_type;

  /// EnzoMethodPpm

  bool                       ppm_diffusion;
  bool                       ppm_dual_energy;
  double                     ppm_dual_energy_eta_1;
  double                     ppm_dual_energy_eta_2;
  int                        ppm_flattening;
  int                        ppm_minimum_pressure_support_parameter;
  double                     ppm_number_density_floor;
  double                     ppm_density_floor;
  double                     ppm_pressure_floor;
  bool                       ppm_pressure_free;
  double                     ppm_temperature_floor;
  bool                       ppm_steepening;
  bool                       ppm_use_minimum_pressure_support;
  double                     ppm_mol_weight;

  double                     field_gamma;

  /// Cosmology
  bool                       physics_cosmology;
  double                     physics_cosmology_hubble_constant_now;
  double                     physics_cosmology_omega_matter_now;
  double                     physics_cosmology_omega_lamda_now;
  double                     physics_cosmology_omega_baryon_now;
  double                     physics_cosmology_omega_cdm_now;
  double                     physics_cosmology_comoving_box_size;
  double                     physics_cosmology_max_expansion_rate;
  double                     physics_cosmology_initial_redshift;
  double                     physics_cosmology_final_redshift;

  /// Gravity
  bool                       physics_gravity;

  /// EnzoInitialCosmology;
  double                     initial_cosmology_temperature;
  
  /// EnzoInitialCollapse
  int                        initial_collapse_rank;
  int                        initial_collapse_array[3];
  double                     initial_collapse_radius_relative;
  double                     initial_collapse_particle_ratio;
  double                     initial_collapse_mass;
  double                     initial_collapse_temperature;

  /// EnzoInitialMusic

  std::vector < std::string > initial_music_field_files;
  std::vector < std::string > initial_music_field_datasets;
  std::vector < std::string > initial_music_field_names;
  std::vector < std::string > initial_music_field_coords;
  
  std::vector < std::string > initial_music_particle_files;
  std::vector < std::string > initial_music_particle_datasets;
  std::vector < std::string > initial_music_particle_coords;
  std::vector < std::string > initial_music_particle_types;
  std::vector < std::string > initial_music_particle_attributes;
  
  /// EnzoInitialPm
  std::string                initial_pm_field;
  double                     initial_pm_mpp;
  int                        initial_pm_level;

  /// EnzoInitialSedovArray[23]
  int                        initial_sedov_rank;
  int                        initial_sedov_array[3];
  double                     initial_sedov_radius_relative;
  double                     initial_sedov_pressure_in;
  double                     initial_sedov_pressure_out;
  double                     initial_sedov_density;

  /// EnzoInitialSedovRandom
  int                        initial_sedov_random_array[3];
  bool                       initial_sedov_random_half_empty; 
  bool                       initial_sedov_random_grackle_cooling;
  int                        initial_sedov_random_max_blasts;
  double                     initial_sedov_random_radius_relative;
  double                     initial_sedov_random_pressure_in;
  double                     initial_sedov_random_pressure_out;
  double                     initial_sedov_random_density;
  int                        initial_sedov_random_te_multiplier;

  /// EnzoInitialSoup
  int                        initial_soup_rank;
  std::string                initial_soup_file;
  bool                       initial_soup_rotate;
  int                        initial_soup_array[3];
  double                     initial_soup_d_pos[3];
  double                     initial_soup_d_size[3];
  double                     initial_soup_pressure_in;
  double                     initial_soup_pressure_out;
  double                     initial_soup_density;

  /// EnzoInitialTurbulence
  double                     initial_turbulence_density;
  double                     initial_turbulence_pressure;
  double                     initial_turbulence_temperature;

  /// EnzoProlong
  std::string                interpolation_method;

  /// EnzoMethodHeat
  double                     method_heat_alpha;

  /// EnzoMethodHydro
  std::string                method_hydro_method;
  bool                       method_hydro_dual_energy;
  double                     method_hydro_dual_energy_eta_1;
  double                     method_hydro_dual_energy_eta_2;
  std::string                method_hydro_reconstruct_method;
  bool                       method_hydro_reconstruct_conservative;
  bool                       method_hydro_reconstruct_positive;
  std::string                method_hydro_riemann_solver;

  /// EnzoMethodNull
  double                     method_null_dt;

  /// EnzoMethodTurbulence
  double                     method_turbulence_edot;
  double                     method_turbulence_mach_number;

  /// EnzoMethodGravity
  double                     method_gravity_grav_const;
  std::string                method_gravity_solver;
  int                        method_gravity_order;
  bool                       method_gravity_accumulate;

  /// EnzoMethodPmDeposit

  double                     method_pm_deposit_alpha;

  /// EnzoMethodPmUpdate

  double                     method_pm_update_max_dt;

  ///==============
  /// EnzoSolverMg0
  ///==============

  /// Solver index for multigrid pre-smoother
  
  std::vector<int>           solver_pre_smooth;

  /// Solver index for multigrid post-smoother

  std::vector<int>           solver_post_smooth;

  /// Solver index for multigrid "last"-smoother

  std::vector<int>           solver_last_smooth;

  /// Solver index for multigrid coarse solver
  
  std::vector<int>           solver_coarse_solve;

  /// Solver index for domain decomposition (dd) domain solver
  
  std::vector<int>           solver_domain_solve;

  /// Weighting factor for smoother
  
  std::vector<double>        solver_weight;

  /// Whether to start the iterative solver using the previous solution

  std::vector<int>           solver_restart_cycle;

  /// EnzoSolver<Krylov>
  
  /// Solver index for Krylov solver preconditioner
  std::vector<int>           solver_precondition;

  /// Whether the solver is for an isolated Block, e.g. for
  /// Mg0 coarse grid solver
  std::vector<int>           solver_local;

  std::vector<int>           solver_coarse_level;
  std::vector<int>           solver_is_unigrid;

  /// Stop at specified redshift for cosmology
  double                     stopping_redshift;

 
#ifdef CONFIG_USE_GRACKLE

  /// EnzoMethodGrackle

  code_units      method_grackle_units;
  chemistry_data  method_grackle_chemistry;

#endif /* CONFIG_USE_GRACKLE */

};

extern EnzoConfig g_enzo_config;

#endif /* PARAMETERS_ENZO_CONFIG_HPP */

