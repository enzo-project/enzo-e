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
      field_uniform_density(1.0),
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
      // FluidProps
      physics_fluid_props_de_config(),
      physics_fluid_props_eos_variant(),
      physics_fluid_props_fluid_floor_config(),
      physics_fluid_props_mol_weight(0.0),


      //--------------------
      // INITIAL [sorted]
      //--------------------

      // EnzoInitialAccretionTest
      initial_accretion_test_gas_density(0.0),
      initial_accretion_test_gas_pressure(0.0),
      initial_accretion_test_gas_radial_velocity(0.0),
      initial_accretion_test_sink_mass(0.0),
      // EnzoInitialBBTest
      initial_bb_test_angular_rotation_velocity(0.0),
      initial_bb_test_external_density(0.0),
      initial_bb_test_fluctuation_amplitude(0.0),
      initial_bb_test_mean_density(0.0),
      initial_bb_test_nominal_sound_speed(0.0),
      initial_bb_test_truncation_radius(0.0),
      // EnzoInitialBCenter
      initial_bcenter_update_etot(false),
      // EnzoInitialBurkertBodenheimer
      initial_burkertbodenheimer_densityprofile(1),
      initial_burkertbodenheimer_mass(0.0),
      initial_burkertbodenheimer_outer_velocity(-1),
      initial_burkertbodenheimer_particle_ratio(0.0),
      initial_burkertbodenheimer_radius_relative(0.0),
      initial_burkertbodenheimer_rank(0),
      initial_burkertbodenheimer_rotating(true),
      initial_burkertbodenheimer_temperature(0.0),
      // EnzoInitialCollapse
      initial_collapse_mass(0.0),
      initial_collapse_particle_ratio(0.0),
      initial_collapse_radius_relative(0.0),
      initial_collapse_rank(0),
      initial_collapse_temperature(0.0),
      // EnzoInitialCosmology
      initial_cosmology_temperature(0.0),
      // EnzoInitialHdf5
      initial_hdf5_blocking(),
      initial_hdf5_field_coords(),
      initial_hdf5_field_levels(),
      initial_hdf5_field_datasets(),
      initial_hdf5_field_files(),
      initial_hdf5_field_names(),
      initial_hdf5_format(),
      initial_hdf5_max_level(),
      initial_hdf5_monitor_iter(),
      initial_hdf5_particle_attributes(),
      initial_hdf5_particle_levels(),
      initial_hdf5_particle_coords(),
      initial_hdf5_particle_datasets(),
      initial_hdf5_particle_files(),
      initial_hdf5_particle_types(),
      //   AE: Maybe these values (and those in cpp) don't matter
      //       are they overwritten by the read-in (even when not found in param file)?
      // EnzoInitialIsolatedGalaxy
      initial_IG_analytic_velocity(false),
      initial_IG_disk_mass(42.9661),            // Gas disk mass in code units
      initial_IG_disk_metal_fraction(1.0E-10),         // Gas disk metal fraction
      initial_IG_disk_temperature(1e4),         // Gas disk temperature in K
      initial_IG_gas_fraction(0.2),             // Gas disk M_gas / M_star
      initial_IG_gas_halo_density(0.0),          // Gas halo uniform density (ignored if zero)
      initial_IG_gas_halo_mass(0.1),             // Gas halo total mass in code units
      initial_IG_gas_halo_metal_fraction(1.0E-10),      // Gas halo metal fraction
      initial_IG_gas_halo_radius(1.0),           // Gas halo maximum radius in code units
      initial_IG_gas_halo_temperature(1e4),      // Gas halo initial temperature
      initial_IG_include_recent_SF(false),
      initial_IG_live_dm_halo(false),
      initial_IG_recent_SF_bin_size(5.0),
      initial_IG_recent_SF_end(0.0),
      initial_IG_recent_SF_seed(12345),
      initial_IG_recent_SF_SFR(2.0),
      initial_IG_recent_SF_start(-100.0),
      initial_IG_scale_height(0.00343218),      // Gas disk scale height in code units
      initial_IG_scale_length(0.0343218),       // Gas disk scale length in code units
      initial_IG_stellar_bulge(false),
      initial_IG_stellar_disk(false),
      initial_IG_use_gas_particles(false),       //
      // EnzoInitialMergeSinksTest
      initial_merge_sinks_test_particle_data_filename(""),
      // EnzoInitialMusic
      initial_music_field_coords(),
      initial_music_field_datasets(),
      initial_music_field_files(),
      initial_music_field_names(),
      initial_music_particle_attributes(),
      initial_music_particle_coords(),
      initial_music_particle_datasets(),
      initial_music_particle_files(),
      initial_music_particle_types(),
      initial_music_throttle_close_count(),
      initial_music_throttle_group_size(),
      initial_music_throttle_internode(),
      initial_music_throttle_intranode(),
      initial_music_throttle_node_files(),
      initial_music_throttle_seconds_delay(),
      initial_music_throttle_seconds_stagger(),
      // EnzoInitialPm
      initial_pm_field(""),
      initial_pm_level(0),
      initial_pm_mpp(0.0),
      initial_sedov_density(0.0),
      // EnzoInitialSedovArray[23]
      initial_sedov_pressure_in(0.0),
      initial_sedov_pressure_out(0.0),
      initial_sedov_radius_relative(0.0),
      // EnzoInitialSedovRandom
      initial_sedov_random_density(0.0),
      initial_sedov_random_grackle_cooling(false),
      initial_sedov_random_half_empty(false),
      initial_sedov_random_max_blasts(0),
      initial_sedov_random_pressure_in(0.0),
      initial_sedov_random_pressure_out(0.0),
      initial_sedov_random_radius_relative(0.0),
      initial_sedov_random_te_multiplier(0),
      initial_sedov_rank(0),
      // EnzoInitialTurbulence
      initial_turbulence_density(0.0),
      initial_turbulence_pressure(0.0),
      initial_turbulence_temperature(0.0),

      //--------------------
      // METHODS [sorted]
      //--------------------

      // EnzoMethodCheck
      method_check_num_files(1),
      method_check_ordering("order_morton"),
      method_check_dir(),
      method_check_monitor_iter(0),
      // EnzoMethodCheckGravity
      method_check_gravity_particle_type(),
      // EnzoMethodGravity
      method_gravity_type_super(),
      // EnzoMethodTurbulence
      method_turbulence_edot(0.0),
      method_turbulence_mach_number(0.0),
      /// EnzoMethodInference
      method_inference_level_base(0),
      method_inference_level_array(0),
      method_inference_level_infer(0),
      method_inference_field_group(),
      method_inference_overdensity_threshold(0),
      // EnzoProlong
      prolong_enzo_type(),
      prolong_enzo_positive(true),
      prolong_enzo_use_linear(false),
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
      solver_coarse_level(),
      solver_is_unigrid(),
      // EnzoStopping
      stopping_redshift()

  {
    for (int axis=0; axis<3; axis++) {
      initial_sedov_array[axis] = 0;
      initial_sedov_random_array[axis] = 0;
      initial_collapse_array[axis] = 0;
      initial_IG_center_position[axis] = 0.5;
      initial_IG_bfield[axis]         = 0.0;
      initial_accretion_test_sink_position[axis] = 0.0;
      initial_accretion_test_sink_velocity[axis] = 0.0;
    }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();

protected: // methods

  void read_adapt_(Parameters *);

  void read_field_(Parameters *);

  //--------------------
  // read_initial [sorted]
  //--------------------
  void read_initial_accretion_test_(Parameters *);
  void read_initial_bb_test_(Parameters *);
  void read_initial_bcenter_(Parameters *);
  void read_initial_burkertbodenheimer_(Parameters *);
  void read_initial_collapse_(Parameters *);
  void read_initial_cosmology_(Parameters *);
  void read_initial_hdf5_(Parameters *);
  void read_initial_isolated_galaxy_(Parameters *);
  void read_initial_merge_sinks_test_(Parameters *);
  void read_initial_music_(Parameters *);
  void read_initial_pm_(Parameters *);
  void read_initial_sedov_(Parameters *);
  void read_initial_sedov_random_(Parameters *);
  void read_initial_turbulence_(Parameters *);

  //--------------------
  // read_method [sorted]
  //--------------------
  void read_method_background_acceleration_(Parameters *);
  void read_method_check_(Parameters *);
  void read_method_gravity_(Parameters *);
  void read_method_inference_(Parameters *);
  void read_method_turbulence_(Parameters *);
  
  void read_physics_(Parameters *);
  void read_physics_fluid_props_(Parameters *);

  void read_prolong_enzo_(Parameters *);

  void read_solvers_(Parameters *);

  void read_stopping_(Parameters *);


public: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Refine

  std::vector <std::string>  adapt_mass_type;

  double                     field_uniform_density;

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

  /// FluidProps
  EnzoDualEnergyConfig       physics_fluid_props_de_config;
  EnzoEOSVariant             physics_fluid_props_eos_variant;
  EnzoFluidFloorConfig       physics_fluid_props_fluid_floor_config;
  double                     physics_fluid_props_mol_weight;

  /// EnzoInitialBCenter;
  bool                       initial_bcenter_update_etot;

  /// EnzoInitialBurkertBodenheimer
  int                        initial_burkertbodenheimer_rank;
  int                        initial_burkertbodenheimer_array[3];
  double                     initial_burkertbodenheimer_radius_relative;
  double                     initial_burkertbodenheimer_particle_ratio;
  double                     initial_burkertbodenheimer_mass;
  double                     initial_burkertbodenheimer_temperature;
  int                        initial_burkertbodenheimer_densityprofile;
  bool                       initial_burkertbodenheimer_rotating;
  double                     initial_burkertbodenheimer_outer_velocity;

  /// EnzoInitialCosmology;
  double                     initial_cosmology_temperature;

  /// EnzoInitialCollapse
  int                        initial_collapse_rank;
  int                        initial_collapse_array[3];
  double                     initial_collapse_radius_relative;
  double                     initial_collapse_particle_ratio;
  double                     initial_collapse_mass;
  double                     initial_collapse_temperature;

  /// EnzoInitialHdf5

  int                         initial_hdf5_max_level;
  std::string                 initial_hdf5_format;
  int                         initial_hdf5_blocking[3];
  int                         initial_hdf5_monitor_iter;
  std::vector < std::string > initial_hdf5_field_files;
  std::vector < std::string > initial_hdf5_field_datasets;
  std::vector < std::string > initial_hdf5_field_names;
  std::vector < std::string > initial_hdf5_field_coords;
  std::vector < int >         initial_hdf5_field_levels;
  std::vector < std::string > initial_hdf5_particle_files;
  std::vector < std::string > initial_hdf5_particle_datasets;
  std::vector < std::string > initial_hdf5_particle_coords;
  std::vector < std::string > initial_hdf5_particle_types;
  std::vector < std::string > initial_hdf5_particle_attributes;
  std::vector < int >         initial_hdf5_particle_levels;

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
  bool                        initial_music_throttle_internode;
  bool                        initial_music_throttle_intranode;
  bool                        initial_music_throttle_node_files;
  int                         initial_music_throttle_close_count;
  int                         initial_music_throttle_group_size;
  double                      initial_music_throttle_seconds_stagger;
  double                      initial_music_throttle_seconds_delay;

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

  /// EnzoInitialTurbulence
  double                     initial_turbulence_density;
  double                     initial_turbulence_pressure;
  double                     initial_turbulence_temperature;

  /// EnzoInitialIsolatedGalaxy
  bool                       initial_IG_analytic_velocity;
  bool                       initial_IG_include_recent_SF;
  bool                       initial_IG_live_dm_halo;
  bool                       initial_IG_stellar_bulge;
  bool                       initial_IG_stellar_disk;
  bool                       initial_IG_use_gas_particles;
  double                     initial_IG_bfield[3];
  double                     initial_IG_center_position[3];
  double                     initial_IG_disk_mass;
  double                     initial_IG_disk_metal_fraction;
  double                     initial_IG_disk_temperature;
  double                     initial_IG_gas_fraction;
  double                     initial_IG_gas_halo_density;
  double                     initial_IG_gas_halo_mass;
  double                     initial_IG_gas_halo_metal_fraction;
  double                     initial_IG_gas_halo_radius;
  double                     initial_IG_gas_halo_temperature;
  double                     initial_IG_recent_SF_bin_size;
  double                     initial_IG_recent_SF_end;
  double                     initial_IG_recent_SF_SFR;
  double                     initial_IG_recent_SF_start;
  double                     initial_IG_scale_height;
  double                     initial_IG_scale_length;
  int                        initial_IG_recent_SF_seed;

  // EnzoInitialMergeSinksTest
  std::string                initial_merge_sinks_test_particle_data_filename;

  // EnzoInitialAccretionTest
  double                     initial_accretion_test_gas_density;
  double                     initial_accretion_test_gas_pressure;
  double                     initial_accretion_test_gas_radial_velocity;
  double                     initial_accretion_test_sink_mass;
  double                     initial_accretion_test_sink_position[3];
  double                     initial_accretion_test_sink_velocity[3];

  // EnzoInitialBBTest
  double                     initial_bb_test_angular_rotation_velocity;
  double                     initial_bb_test_center[3];
  double                     initial_bb_test_drift_velocity[3];
  double                     initial_bb_test_external_density;
  double                     initial_bb_test_fluctuation_amplitude;
  double                     initial_bb_test_mean_density;
  double                     initial_bb_test_nominal_sound_speed;
  double                     initial_bb_test_truncation_radius;

  //--------------------
  // EnzoMethod
  //--------------------

  /// EnzoMethodCheck
  int                        method_check_num_files;
  std::string                method_check_ordering;
  std::vector<std::string>   method_check_dir;
  int                        method_check_monitor_iter;
  bool                       method_check_include_ghosts;

  /// EnzoMethodCheckGravity
  std::string                method_check_gravity_particle_type;

  /// EnzoMethodGravity
  std::string                method_gravity_type_super;

  /// EnzoMethodInference
  int                        method_inference_level_base;
  int                        method_inference_level_array;
  int                        method_inference_level_infer;
  std::string                method_inference_field_group;
  float                      method_inference_overdensity_threshold;

  /// EnzoMethodTurbulence
  double                     method_turbulence_edot;
  double                     method_turbulence_mach_number;

  std::string                prolong_enzo_type;
  bool                       prolong_enzo_positive;
  bool                       prolong_enzo_use_linear;

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

  /// Mg0 coarse grid solver

  std::vector<int>           solver_coarse_level;
  std::vector<int>           solver_is_unigrid;

  /// Stop at specified redshift for cosmology
  double                     stopping_redshift;

};

extern EnzoConfig g_enzo_config;

#endif /* PARAMETERS_ENZO_CONFIG_HPP */
