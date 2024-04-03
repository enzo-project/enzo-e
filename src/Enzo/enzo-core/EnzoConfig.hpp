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
      ppm_flattening(0),
      ppm_minimum_pressure_support_parameter(0),
      ppm_pressure_free(false),
      ppm_steepening(false),
      ppm_use_minimum_pressure_support(false),
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
      // Gravity
      physics_gravity_grav_constant_codeU(-1.0),

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
      // EnzoInitialCloud
      initial_cloud_center_x(0.0),
      initial_cloud_center_y(0.0),
      initial_cloud_center_z(0.0),
      initial_cloud_density_cloud(0.0),
      initial_cloud_density_wind(0.0),
      initial_cloud_eint_wind(0.0),
      initial_cloud_etot_wind(0.0),
      initial_cloud_initialize_uniform_bfield(false),
      initial_cloud_metal_mass_frac(0.0),
      initial_cloud_perturb_Nwaves(0),
      initial_cloud_perturb_amplitude(0.0),
      initial_cloud_perturb_min_wavelength(std::numeric_limits<double>::min()),
      initial_cloud_perturb_max_wavelength(std::numeric_limits<double>::min()),
      initial_cloud_perturb_seed(0),
      initial_cloud_radius(0.),
      initial_cloud_subsample_n(0),
      initial_cloud_velocity_wind(0.0),
      // EnzoInitialCollapse
      initial_collapse_mass(0.0),
      initial_collapse_particle_ratio(0.0),
      initial_collapse_radius_relative(0.0),
      initial_collapse_rank(0),
      initial_collapse_temperature(0.0),
      // EnzoInitialCosmology
      initial_cosmology_temperature(0.0),
      // EnzoInitialFeedbackTest
      initial_feedback_test_density(),
      initial_feedback_test_e_density(),
      initial_feedback_test_from_file(),
      initial_feedback_test_HeI_density(),
      initial_feedback_test_HeII_density(),
      initial_feedback_test_HeIII_density(),
      initial_feedback_test_HI_density(),
      initial_feedback_test_HII_density(),
      initial_feedback_test_metal_fraction(),
      initial_feedback_test_star_mass(),
      initial_feedback_test_temperature(),
      initial_feedback_test_luminosity(),
      // EnzoGrackleTest
      initial_grackle_test_maximum_H_number_density(1000.0),
      initial_grackle_test_maximum_metallicity(1.0),
      initial_grackle_test_maximum_temperature(1.0E8),
      initial_grackle_test_minimum_H_number_density(0.1),
      initial_grackle_test_minimum_metallicity(1.0E-4),
      initial_grackle_test_minimum_temperature(10.0),
      // EnzoInitialHdf5
      initial_hdf5_blocking(),
      initial_hdf5_field_coords(),
      initial_hdf5_field_datasets(),
      initial_hdf5_field_files(),
      initial_hdf5_field_names(),
      initial_hdf5_format(),
      initial_hdf5_max_level(),
      initial_hdf5_monitor_iter(),
      initial_hdf5_particle_attributes(),
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
      // EnzoInitialShuCollapse
      initial_shu_collapse_central_sink_exists(false),
      initial_shu_collapse_central_sink_mass(0.0),
      initial_shu_collapse_external_density(0.0),
      initial_shu_collapse_instability_parameter(0.0),
      initial_shu_collapse_nominal_sound_speed(0.0),
      initial_shu_collapse_truncation_radius(0.0),
      // EnzoInitialSoup
      initial_soup_density(0.0),
      initial_soup_file(""),
      initial_soup_pressure_in(0.0),
      initial_soup_pressure_out(0.0),
      initial_soup_rank(0),
      initial_soup_rotate(false),
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
      method_check_include_ghosts(false),
      /// EnzoMethodFeedback
      method_feedback_ejecta_mass(0.0),
      method_feedback_ejecta_metal_fraction(0.0),
      method_feedback_flavor(""),
      method_feedback_ke_fraction(0.0),
      method_feedback_radius(-1.0),
      method_feedback_shift_cell_center(true),
      method_feedback_stencil(3),
      method_feedback_supernova_energy(1.0),
      method_feedback_time_first_sn(-1.0), // in Myr
      method_feedback_use_ionization_feedback(false),
      /// EnzoMethodFeedbackSTARSS
      method_feedback_supernovae(true),
      method_feedback_unrestricted_sn(true),
      method_feedback_stellar_winds(true),
      method_feedback_min_level(0),
      method_feedback_analytic_SNR_shell_mass(true),
      method_feedback_fade_SNR(true),
      method_feedback_NEvents(-1),
      // EnzoMethodCheckGravity
      method_check_gravity_particle_type(),
      /// EnzoMethodStarMaker
      method_star_maker_flavor(""),
      method_star_maker_use_density_threshold(false),           // check above density threshold before SF
      method_star_maker_use_velocity_divergence(false),         // check for converging flow before SF
      method_star_maker_use_dynamical_time(false),              //
      method_star_maker_use_altAlpha(false), // alternate virial parameter calculation
      method_star_maker_use_cooling_time(false), 
      method_star_maker_use_self_gravitating(false),           //
      method_star_maker_use_h2_self_shielding(false),
      method_star_maker_use_jeans_mass(false),
      method_star_maker_use_overdensity_threshold(false),
      method_star_maker_use_critical_metallicity(false),
      method_star_maker_use_temperature_threshold(false),
      method_star_maker_critical_metallicity(0.0),
      method_star_maker_temperature_threshold(1.0E4),
      method_star_maker_number_density_threshold(0.0),      // Number density threshold in cgs
      method_star_maker_maximum_mass_fraction(0.05),            // maximum cell mass fraction to convert to stars
      method_star_maker_efficiency(0.01),            // star maker efficiency
      method_star_maker_minimum_star_mass(0.0),    // minium star particle mass in solar masses
      method_star_maker_maximum_star_mass(-1.0),    // maximum star particle mass in solar masses
      method_star_maker_min_level(0), // minimum refinement level for star formation
      method_star_maker_turn_off_probability(false),
      // EnzoMethodM1Closure
      method_m1_closure(false),
      method_m1_closure_N_groups(1), // # of frequency bins
      method_m1_closure_flux_function("GLF"), // which flux function to use
      method_m1_closure_hll_file("hll_evals.list"), // path to hll eigenvalue table
      method_m1_closure_clight_frac(1.0), // fraction of speed of light value to use
      method_m1_closure_photon_escape_fraction(1.0),
      method_m1_closure_radiation_spectrum("custom"), // Type of radiation spectrum to use for star particles 
      method_m1_closure_temperature_blackbody(-1.0), // requires radiation_spectrum="blackbody"
      method_m1_closure_particle_luminosity(-1.0), // specify emission rate from sources
      method_m1_closure_SED(), // supply list of emission rates for all groups (radiation_spectrum="custom")
      method_m1_closure_min_photon_density(0.0),
      method_m1_closure_attenuation(true),
      method_m1_closure_thermochemistry(true),
      method_m1_closure_recombination_radiation(false),
      method_m1_closure_H2_photodissociation(false),
      method_m1_closure_lyman_werner_background(false),
      method_m1_closure_LWB_J21(-1.0),
      method_m1_closure_cross_section_calculator("vernier"), // what type of cross section calculator to use ("vernier", "vernier_average", "custom")
      method_m1_closure_sigmaN(), // user-defined cross sections (requires cross_section_calculator = "custom")
      method_m1_closure_sigmaE(), 
      method_m1_closure_energy_lower(),
      method_m1_closure_energy_upper(),
      method_m1_closure_energy_mean(), // mean energy for bin
      // EnzoMethodTurbulence
      method_turbulence_edot(0.0),
      method_turbulence_mach_number(0.0),
      // EnzoMethodGrackle
      method_grackle_use_grackle(false),
      method_grackle_chemistry(),
      method_grackle_use_cooling_timestep(false),
      method_grackle_radiation_redshift(-1.0),
      // EnzoMethodGravity
      method_gravity_solver(""),
      method_gravity_order(4),
      method_gravity_dt_max(1.0e10),
      method_gravity_accumulate(false),
      // EnzoMethodBackgroundAcceleration
      method_background_acceleration_flavor(""),
      method_background_acceleration_mass(0.0),
      method_background_acceleration_DM_mass(0.0),
      method_background_acceleration_bulge_mass(0.0),
      method_background_acceleration_core_radius(0.0),
      method_background_acceleration_bulge_radius(0.0),
      method_background_acceleration_stellar_mass(0.0),
      method_background_acceleration_DM_mass_radius(0.0),
      method_background_acceleration_stellar_scale_height_r(0.0),
      method_background_acceleration_stellar_scale_height_z(0.0),
      method_background_acceleration_apply_acceleration(true),
      // EnzoMethodMHDVlct
      method_vlct_riemann_solver(""),
      method_vlct_time_scheme(""),
      method_vlct_reconstruct_method(""),
      method_vlct_theta_limiter(0.0),
      method_vlct_mhd_choice(""),
      // EnzoMethodMergeSinks
      method_merge_sinks_merging_radius_cells(0.0),
      // EnzoMethodAccretion
      method_accretion_accretion_radius_cells(0.0),
      method_accretion_flavor(""),
      method_accretion_physical_density_threshold_cgs(0.0),
      method_accretion_max_mass_fraction(0.0),
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
      initial_cloud_uniform_bfield[axis] = 0;
      initial_sedov_array[axis] = 0;
      initial_sedov_random_array[axis] = 0;
      initial_soup_array[axis] = 0;
      initial_soup_d_pos[axis] = 0;
      initial_soup_d_size[axis] = 0;
      initial_collapse_array[axis] = 0;
      initial_IG_center_position[axis] = 0.5;
      initial_IG_bfield[axis]         = 0.0;
      initial_accretion_test_sink_position[axis] = 0.0;
      initial_accretion_test_sink_velocity[axis] = 0.0;
      method_background_acceleration_center[axis] = 0.5;
      method_background_acceleration_angular_momentum[axis] = 0;

      initial_feedback_test_position[axis] = 0.5;
    }
    method_background_acceleration_angular_momentum[2] = 1;
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
  void read_initial_cloud_(Parameters *);
  void read_initial_collapse_(Parameters *);
  void read_initial_cosmology_(Parameters *);
  void read_initial_feedback_test_(Parameters *);
  void read_initial_grackle_(Parameters *);
  void read_initial_hdf5_(Parameters *);
  void read_initial_isolated_galaxy_(Parameters *);
  void read_initial_merge_sinks_test_(Parameters *);
  void read_initial_music_(Parameters *);
  void read_initial_pm_(Parameters *);
  void read_initial_sedov_(Parameters *);
  void read_initial_sedov_random_(Parameters *);
  void read_initial_shu_collapse_(Parameters *);
  void read_initial_soup_(Parameters *);
  void read_initial_turbulence_(Parameters *);

  //--------------------
  // read_method [sorted]
  //--------------------
  void read_method_accretion_(Parameters *);
  void read_method_background_acceleration_(Parameters *);
  void read_method_check_(Parameters *);
  void read_method_feedback_(Parameters *);
  void read_method_grackle_(Parameters *);
  void read_method_gravity_(Parameters *);
  void read_method_merge_sinks_(Parameters *);
  void read_method_ppm_(Parameters *);
  void read_method_star_maker_(Parameters *);
  void read_method_m1_closure_(Parameters *);
  void read_method_turbulence_(Parameters *);
  void read_method_vlct_(Parameters *);
  
  void read_physics_(Parameters *);
  void read_physics_fluid_props_(Parameters *);
  void read_physics_gravity_(Parameters *);

  void read_prolong_enzo_(Parameters *);

  void read_solvers_(Parameters *);

  void read_stopping_(Parameters *);


public: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Refine

  std::vector <std::string>  adapt_mass_type;

  /// EnzoMethodPpm

  bool                       ppm_diffusion;
  int                        ppm_flattening;
  int                        ppm_minimum_pressure_support_parameter;
  bool                       ppm_pressure_free;
  bool                       ppm_steepening;
  bool                       ppm_use_minimum_pressure_support;

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

  /// Gravity
  double                     physics_gravity_grav_constant_codeU;

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

  /// EnzoInitialCloud;
  double                     initial_cloud_center_x;
  double                     initial_cloud_center_y;
  double                     initial_cloud_center_z;
  double                     initial_cloud_density_cloud;
  double                     initial_cloud_density_wind;
  double                     initial_cloud_eint_wind;
  double                     initial_cloud_etot_wind;
  bool                       initial_cloud_initialize_uniform_bfield;
  double                     initial_cloud_uniform_bfield[3];
  double                     initial_cloud_metal_mass_frac;
  int                        initial_cloud_perturb_Nwaves;
  double                     initial_cloud_perturb_amplitude;
  double                     initial_cloud_perturb_min_wavelength;
  double                     initial_cloud_perturb_max_wavelength;
  unsigned int               initial_cloud_perturb_seed;
  double                     initial_cloud_radius;
  int                        initial_cloud_subsample_n;
  double                     initial_cloud_velocity_wind;

  /// EnzoInitialCosmology;
  double                     initial_cosmology_temperature;

  /// EnzoInitialCollapse
  int                        initial_collapse_rank;
  int                        initial_collapse_array[3];
  double                     initial_collapse_radius_relative;
  double                     initial_collapse_particle_ratio;
  double                     initial_collapse_mass;
  double                     initial_collapse_temperature;

  /// EnzoGrackleTest
  double                     initial_grackle_test_maximum_H_number_density;
  double                     initial_grackle_test_maximum_metallicity;
  double                     initial_grackle_test_maximum_temperature;
  double                     initial_grackle_test_minimum_H_number_density;
  double                     initial_grackle_test_minimum_metallicity;
  double                     initial_grackle_test_minimum_temperature;

  /// EnzoInitialHdf5

  int                         initial_hdf5_max_level;
  std::string                 initial_hdf5_format;
  int                         initial_hdf5_blocking[3];
  int                         initial_hdf5_monitor_iter;
  std::vector < std::string > initial_hdf5_field_files;
  std::vector < std::string > initial_hdf5_field_datasets;
  std::vector < std::string > initial_hdf5_field_names;
  std::vector < std::string > initial_hdf5_field_coords;
  std::vector < std::string > initial_hdf5_particle_files;
  std::vector < std::string > initial_hdf5_particle_datasets;
  std::vector < std::string > initial_hdf5_particle_coords;
  std::vector < std::string > initial_hdf5_particle_types;
  std::vector < std::string > initial_hdf5_particle_attributes;

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

  /// EnzoInitialFeedbackTest

  double                     initial_feedback_test_position[3];
  double                     initial_feedback_test_luminosity;
  double                     initial_feedback_test_density;
  double                     initial_feedback_test_HI_density;
  double                     initial_feedback_test_HII_density;
  double                     initial_feedback_test_HeI_density;
  double                     initial_feedback_test_HeII_density;
  double                     initial_feedback_test_HeIII_density;
  double                     initial_feedback_test_e_density;
  double                     initial_feedback_test_star_mass;
  double                     initial_feedback_test_temperature;
  bool                       initial_feedback_test_from_file;
  double                     initial_feedback_test_metal_fraction;

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

  // EnzoInitialShuCollapse
  bool                       initial_shu_collapse_central_sink_exists;
  double                     initial_shu_collapse_center[3];
  double                     initial_shu_collapse_central_sink_mass;
  double                     initial_shu_collapse_drift_velocity[3];
  double                     initial_shu_collapse_external_density;
  double                     initial_shu_collapse_instability_parameter;
  double                     initial_shu_collapse_nominal_sound_speed;
  double                     initial_shu_collapse_truncation_radius;

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

  /// EnzoMethodFeedback
  std::string               method_feedback_flavor;
  double                    method_feedback_ejecta_mass;
  double                    method_feedback_supernova_energy;
  double                    method_feedback_ejecta_metal_fraction;
  double                    method_feedback_ke_fraction;
  double                    method_feedback_time_first_sn;
  int                       method_feedback_stencil;
  double                    method_feedback_radius;
  bool                      method_feedback_shift_cell_center;
  bool                      method_feedback_use_ionization_feedback;

  
  /// EnzoMethodFeedbackSTARSS
  bool                       method_feedback_supernovae;
  bool                       method_feedback_unrestricted_sn;
  bool                       method_feedback_stellar_winds;
  bool                       method_feedback_radiation;
  int                        method_feedback_min_level;
  bool                       method_feedback_analytic_SNR_shell_mass;
  bool                       method_feedback_fade_SNR;
  int                        method_feedback_NEvents;

  /// EnzoMethodStarMaker
  std::string               method_star_maker_flavor;
  bool                      method_star_maker_use_altAlpha;
  bool                      method_star_maker_use_density_threshold;
  bool                      method_star_maker_use_overdensity_threshold;
  bool                      method_star_maker_use_temperature_threshold;
  bool                      method_star_maker_use_critical_metallicity;
  bool                      method_star_maker_use_velocity_divergence;
  bool                      method_star_maker_use_cooling_time;
  bool                      method_star_maker_use_dynamical_time;
  bool                      method_star_maker_use_h2_self_shielding;
  bool                      method_star_maker_use_jeans_mass;
  bool                      method_star_maker_use_self_gravitating;
  double                    method_star_maker_number_density_threshold;
  double                    method_star_maker_overdensity_threshold;
  double                    method_star_maker_temperature_threshold;
  double                    method_star_maker_critical_metallicity;
  double                    method_star_maker_maximum_mass_fraction;
  double                    method_star_maker_efficiency;
  double                    method_star_maker_minimum_star_mass;
  double                    method_star_maker_maximum_star_mass;
  int                       method_star_maker_min_level;
  bool                      method_star_maker_turn_off_probability;


  /// EnzoMethodM1Closure
  
  bool                      method_m1_closure; 
  int                       method_m1_closure_N_groups;
  std::string               method_m1_closure_flux_function;
  std::string               method_m1_closure_hll_file;
  double                    method_m1_closure_clight_frac;
  double                    method_m1_closure_photon_escape_fraction;
  std::string               method_m1_closure_radiation_spectrum;
  double                    method_m1_closure_temperature_blackbody; 
  double                    method_m1_closure_particle_luminosity;
  std::vector<double>       method_m1_closure_SED;
  double                    method_m1_closure_min_photon_density;
  bool                      method_m1_closure_attenuation;
  bool                      method_m1_closure_thermochemistry;
  bool                      method_m1_closure_recombination_radiation;
  bool                      method_m1_closure_H2_photodissociation;
  bool                      method_m1_closure_lyman_werner_background;
  double                    method_m1_closure_LWB_J21;
  std::string               method_m1_closure_cross_section_calculator;
  std::vector<double>       method_m1_closure_sigmaN;
  std::vector<double>       method_m1_closure_sigmaE;
  std::vector<double>       method_m1_closure_energy_lower;
  std::vector<double>       method_m1_closure_energy_upper;
  std::vector<double>       method_m1_closure_energy_mean;

  /// EnzoMethodTurbulence
  double                     method_turbulence_edot;
  double                     method_turbulence_mach_number;

  /// EnzoMethodGrackle
  bool                       method_grackle_use_grackle;
  GrackleChemistryData       method_grackle_chemistry;
  bool                       method_grackle_use_cooling_timestep;
  double                     method_grackle_radiation_redshift;

  /// EnzoMethodGravity
  std::string                method_gravity_solver;
  int                        method_gravity_order;
  double                     method_gravity_dt_max;
  bool                       method_gravity_accumulate;

  /// EnzoMethodBackgroundAcceleration

  std::string                method_background_acceleration_flavor;
  double                     method_background_acceleration_mass;
  double                     method_background_acceleration_DM_mass;
  double                     method_background_acceleration_bulge_mass;
  double                     method_background_acceleration_core_radius;
  double                     method_background_acceleration_bulge_radius;
  double                     method_background_acceleration_stellar_mass;
  double                     method_background_acceleration_DM_mass_radius;
  double                     method_background_acceleration_stellar_scale_height_r;
  double                     method_background_acceleration_stellar_scale_height_z;
  double                     method_background_acceleration_center[3];
  double                     method_background_acceleration_angular_momentum[3];
  bool                       method_background_acceleration_apply_acceleration;

  /// EnzoMethodMHDVlct
  std::string                method_vlct_riemann_solver;
  std::string                method_vlct_time_scheme;
  std::string                method_vlct_reconstruct_method;
  double                     method_vlct_theta_limiter;
  std::string                method_vlct_mhd_choice;

  /// EnzoMethodMergeSinks
  double                     method_merge_sinks_merging_radius_cells;

  /// EnzoMethodAccretion
  double                     method_accretion_accretion_radius_cells;
  std::string                method_accretion_flavor;
  double                     method_accretion_physical_density_threshold_cgs;
  double                     method_accretion_max_mass_fraction;
  
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
