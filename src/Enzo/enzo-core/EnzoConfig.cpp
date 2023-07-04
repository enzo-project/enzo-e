// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoConfig.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-03
/// @brief    Implementation of the EnzoConfig class

#include "cello.hpp"
#include "enzo.hpp"

extern CProxy_EnzoSimulation proxy_enzo_simulation;

//----------------------------------------------------------------------

EnzoConfig g_enzo_config;

EnzoConfig::EnzoConfig() throw ()
  :
  adapt_mass_type(0),
  ppm_diffusion(false),
  ppm_flattening(0),
  ppm_minimum_pressure_support_parameter(0),
  ppm_pressure_free(false),
  ppm_steepening(false),
  ppm_use_minimum_pressure_support(false),
  field_uniform_density(1.0),
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
  physics_gravity(false),
  // EnzoInitialBCenter
  initial_bcenter_update_etot(false),
  // EnzoInitialBurkertBodenheimer
  initial_burkertbodenheimer_rank(0),
  initial_burkertbodenheimer_radius_relative(0.0),
  initial_burkertbodenheimer_particle_ratio(0.0),
  initial_burkertbodenheimer_mass(0.0),
  initial_burkertbodenheimer_temperature(0.0),
  initial_burkertbodenheimer_densityprofile(1),
  initial_burkertbodenheimer_rotating(true),
  initial_burkertbodenheimer_outer_velocity(-1),
  // EnzoInitialCloud
  initial_cloud_subsample_n(0),
  initial_cloud_radius(0.),
  initial_cloud_center_x(0.0),
  initial_cloud_center_y(0.0),
  initial_cloud_center_z(0.0),
  initial_cloud_density_cloud(0.0),
  initial_cloud_density_wind(0.0),
  initial_cloud_velocity_wind(0.0),
  initial_cloud_etot_wind(0.0),
  initial_cloud_eint_wind(0.0),
  initial_cloud_metal_mass_frac(0.0),
  initial_cloud_initialize_uniform_bfield(false),
  initial_cloud_perturb_stddev(0.0),
  initial_cloud_trunc_dev(0.0),
  initial_cloud_perturb_seed(0),
  // EnzoInitialCosmology
  initial_cosmology_temperature(0.0),
  // EnzoInitialCollapse
  initial_collapse_rank(0),
  initial_collapse_radius_relative(0.0),
  initial_collapse_particle_ratio(0.0),
  initial_collapse_mass(0.0),
  initial_collapse_temperature(0.0),
  // EnzoInitialFeedbackTest
  initial_feedback_test_density(),
  initial_feedback_test_HI_density(),
  initial_feedback_test_HII_density(),
  initial_feedback_test_HeI_density(),
  initial_feedback_test_HeII_density(),
  initial_feedback_test_HeIII_density(),
  initial_feedback_test_e_density(),
  initial_feedback_test_star_mass(),
  initial_feedback_test_temperature(),
  initial_feedback_test_from_file(),
  initial_feedback_test_metal_fraction(0.01),
  initial_feedback_test_luminosity(),
  // EnzoInitialGrackleTest
  initial_grackle_test_maximum_H_number_density(1000.0),
  initial_grackle_test_maximum_metallicity(1.0),
  initial_grackle_test_maximum_temperature(1.0E8),
  initial_grackle_test_minimum_H_number_density(0.1),
  initial_grackle_test_minimum_metallicity(1.0E-4),
  initial_grackle_test_minimum_temperature(10.0),
  initial_grackle_test_reset_energies(0),
  // EnzoInitialHdf5
  initial_hdf5_max_level(),
  initial_hdf5_format(),
  initial_hdf5_blocking(),
  initial_hdf5_monitor_iter(),
  initial_hdf5_field_files(),
  initial_hdf5_field_datasets(),
  initial_hdf5_field_names(),
  initial_hdf5_field_coords(),
  initial_hdf5_particle_files(),
  initial_hdf5_particle_datasets(),
  initial_hdf5_particle_coords(),
  initial_hdf5_particle_types(),
  initial_hdf5_particle_attributes(),
  // EnzoInitialInclinedWave
  initial_inclinedwave_alpha(0.0),
  initial_inclinedwave_beta(0.0),
  initial_inclinedwave_amplitude(0.0),
  initial_inclinedwave_lambda(0.0),
  initial_inclinedwave_parallel_vel(std::numeric_limits<double>::min()),
  initial_inclinedwave_positive_vel(true),
  initial_inclinedwave_wave_type(""),
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
  initial_music_throttle_internode(),
  initial_music_throttle_intranode(),
  initial_music_throttle_node_files(),
  initial_music_throttle_close_count(),
  initial_music_throttle_group_size(),
  initial_music_throttle_seconds_stagger(),
  initial_music_throttle_seconds_delay(),
  // EnzoInitialPm
  initial_pm_field(""),
  initial_pm_mpp(0.0),
  initial_pm_level(0),
  // EnzoInitialSedov[23]
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
  // EnzoInitialShockTube
  initial_shock_tube_setup_name(""),
  initial_shock_tube_aligned_ax(""),
  initial_shock_tube_axis_velocity(0.0),
  initial_shock_tube_trans_velocity(0.0),
  initial_shock_tube_flip_initialize(false),
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
  initial_IG_use_gas_particles(false),      // Set up gas by depositing baryonic particles to grid
  // EnzoMethodCheck
  method_check_num_files(1),
  method_check_ordering("order_morton"),
  method_check_dir(),
  method_check_monitor_iter(0),
  // EnzoInitialMergeSinksTest
  initial_merge_sinks_test_particle_data_filename(""),
  // EnzoInitialAccretionTest
  initial_accretion_test_sink_mass(0.0),
  initial_accretion_test_gas_density(0.0),
  initial_accretion_test_gas_pressure(0.0),
  initial_accretion_test_gas_radial_velocity(0.0),
  // EnzoInitialShuCollapse
  initial_shu_collapse_truncation_radius(0.0),
  initial_shu_collapse_nominal_sound_speed(0.0),
  initial_shu_collapse_instability_parameter(0.0),
  initial_shu_collapse_external_density(0.0),
  initial_shu_collapse_central_sink_exists(false),
  initial_shu_collapse_central_sink_mass(0.0),
  // EnzoInitialBBTest
  initial_bb_test_mean_density(0.0),
  initial_bb_test_fluctuation_amplitude(0.0),
  initial_bb_test_truncation_radius(0.0),
  initial_bb_test_nominal_sound_speed(0.0),
  initial_bb_test_angular_rotation_velocity(0.0),
  initial_bb_test_external_density(0.0),
  // EnzoMethodHeat
  method_heat_alpha(0.0),
  // EnzoMethodHydro
  method_hydro_method(""),
  method_hydro_dual_energy(false),
  method_hydro_dual_energy_eta_1(0.0),
  method_hydro_dual_energy_eta_2(0.0),
  method_hydro_reconstruct_method(""),
  method_hydro_reconstruct_conservative(0),
  method_hydro_reconstruct_positive(0),
  method_hydro_riemann_solver(""),
  // EnzoMethodFeedback,
  method_feedback_flavor(""),
  method_feedback_ejecta_mass(0.0),
  method_feedback_supernova_energy(1.0),
  method_feedback_ejecta_metal_fraction(0.0),
  method_feedback_stencil(3),
  method_feedback_radius(-1),
  method_feedback_shift_cell_center(true),
  method_feedback_ke_fraction(0.0),
  method_feedback_use_ionization_feedback(false),
  method_feedback_time_first_sn(-1), // in Myr
  // EnzoMethodFeedbackSTARSS,
  method_feedback_supernovae(true),
  method_feedback_unrestricted_sn(true),
  method_feedback_stellar_winds(true),
  method_feedback_min_level(0),
  method_feedback_analytic_SNR_shell_mass(true),
  method_feedback_fade_SNR(true),
  method_feedback_NEvents(-1),
  method_feedback_radiation(true),
  // EnzoMethodM1Closure
  method_m1_closure(false),
  method_m1_closure_N_groups(1), // # of frequency bins
  method_m1_closure_flux_function("GLF"), // which flux function to use
  method_m1_closure_hll_file("hll_evals.list"),
  method_m1_closure_clight_frac(1.0), // reduced speed of light value to use
  method_m1_closure_photon_escape_fraction(1.0),
  method_m1_closure_radiation_spectrum("custom"), // Type of radiation spectrum to use for star particles
  method_m1_closure_temperature_blackbody(0.0),
  method_m1_closure_particle_luminosity(-1.0), // Set emission rate for star particles
  method_m1_closure_SED(), // supply list of emission rate fraction for all groups
  method_m1_closure_min_photon_density(0.0),
  method_m1_closure_attenuation(true),
  method_m1_closure_thermochemistry(true),
  method_m1_closure_recombination_radiation(false),
  method_m1_closure_H2_photodissociation(false),
  method_m1_closure_lyman_werner_background(false),
  method_m1_closure_LWB_J21(-1.0),
  method_m1_closure_cross_section_calculator("vernier"),
  method_m1_closure_sigmaN(),
  method_m1_closure_sigmaE(),
  method_m1_closure_energy_lower(),
  method_m1_closure_energy_upper(),
  method_m1_closure_energy_mean(),
  // EnzoMethodStarMaker,
  method_star_maker_flavor(""),                              // star maker type to use
  method_star_maker_use_altAlpha(false),
  method_star_maker_use_density_threshold(false),           // check above density threshold before SF
  method_star_maker_use_velocity_divergence(false),         // check for converging flow before SF
  method_star_maker_use_dynamical_time(false),              // compute t_ff / t_dyn. Otherwise take as 1.0
  method_star_maker_use_cooling_time(false),                // check if t_cool < t_dyn
  method_star_maker_use_overdensity_threshold(false),
  method_star_maker_use_temperature_threshold(false),
  method_star_maker_use_self_gravitating(false),            //
  method_star_maker_use_h2_self_shielding(false),
  method_star_maker_use_jeans_mass(false),
  method_star_maker_number_density_threshold(0.0),         // Number density threshold in cgs
  method_star_maker_overdensity_threshold(0.0),
  method_star_maker_critical_metallicity(0.0),
  method_star_maker_temperature_threshold(1.0E4),
  method_star_maker_maximum_mass_fraction(0.05),            // maximum cell mass fraction to convert to stars
  method_star_maker_efficiency(0.01),            // star maker efficiency per free fall time
  method_star_maker_minimum_star_mass(0.0),    // minimum star particle mass in solar masses
  method_star_maker_maximum_star_mass(-1.0),    // maximum star particle mass in solar masses
  method_star_maker_min_level(0), // minimum AMR level for star formation
  method_star_maker_turn_off_probability(false),
  // EnzoMethodTurbulence
  method_turbulence_edot(0.0),
  method_turbulence_mach_number(0.0),
  method_grackle_use_grackle(false),
  method_grackle_chemistry(),
  method_grackle_use_cooling_timestep(false),
  method_grackle_radiation_redshift(-1.0),
  // EnzoMethodGravity
  method_gravity_grav_const(0.0),
  method_gravity_solver(""),
  method_gravity_order(4),
  method_gravity_dt_max(0.0),
  method_gravity_accumulate(false),
  /// EnzoMethodBackgroundAcceleration
  method_background_acceleration_flavor(""),
  method_background_acceleration_mass(0.0),
  method_background_acceleration_DM_mass(0.0),
  method_background_acceleration_bulge_mass(0.0),
  method_background_acceleration_core_radius(1.0E-10),
  method_background_acceleration_bulge_radius(1.0E-10),
  method_background_acceleration_stellar_mass(0.0),
  method_background_acceleration_DM_mass_radius(0.0),
  method_background_acceleration_stellar_scale_height_r(1.0E-10),
  method_background_acceleration_stellar_scale_height_z(1.0E-10),
  method_background_acceleration_apply_acceleration(true), // for debugging
  /// EnzoMethodPmDeposit
  method_pm_deposit_alpha(0.5),
  /// EnzoMethodPmUpdate
  method_pm_update_max_dt(std::numeric_limits<double>::max()),
  /// EnzoMethodMHDVlct
  method_vlct_riemann_solver(""),
  method_vlct_half_dt_reconstruct_method(""),
  method_vlct_full_dt_reconstruct_method(""),
  method_vlct_theta_limiter(0.0),
  method_vlct_mhd_choice(""),
  /// EnzoMethodMergeSinks
  method_merge_sinks_merging_radius_cells(0.0),
  /// EnzoMethodAccretion
  method_accretion_accretion_radius_cells(0.0),
  method_accretion_flavor(""),
  method_accretion_physical_density_threshold_cgs(0.0),
  method_accretion_max_mass_fraction(0.0),
  /// EnzoMethodSinkMaker
  method_sink_maker_jeans_length_resolution_cells(0.0),
  method_sink_maker_physical_density_threshold_cgs(0.0),
  method_sink_maker_check_density_maximum(false),
  method_sink_maker_max_mass_fraction(0.0),
  method_sink_maker_min_sink_mass_solar(0.0),
  method_sink_maker_max_offset_cell_fraction(0.0),
  method_sink_maker_offset_seed_shift(0),
  /// EnzoProlong
  prolong_enzo_type(),
  prolong_enzo_positive(true),
  prolong_enzo_use_linear(false),
  /// EnzoSolverMg0
  solver_pre_smooth(),
  solver_post_smooth(),
  solver_last_smooth(),
  solver_coarse_solve(),
  solver_domain_solve(),
  solver_weight(),
  solver_restart_cycle(),
  /// EnzoSolver<Krylov>
  solver_precondition(),
  solver_coarse_level(),
  solver_is_unigrid(),
  stopping_redshift()

{
  for (int i=0; i<3; i++) {
    initial_cloud_uniform_bfield[i] = 0;
    initial_sedov_array[i] = 0;
    initial_soup_array[i]  = 0;
    initial_soup_d_pos[i]  = 0.0;
    initial_soup_d_size[i] = 0.0;
    initial_collapse_array[i] = 0;
    initial_IG_center_position[i] = 0.5;
    initial_IG_bfield[i] = 0.0;
    method_background_acceleration_center[i] = 0.5;
    method_background_acceleration_angular_momentum[i] = 0;

    initial_feedback_test_position[i] = 0.5;
  }

  method_background_acceleration_angular_momentum[2] = 1;
}

//----------------------------------------------------------------------

EnzoConfig::~EnzoConfig() throw ()
{ }

//----------------------------------------------------------------------

void EnzoConfig::pup (PUP::er &p)
{

  Config::pup(p);
  TRACEPUP;

  // NOTE: change this function whenever attributes change

  p | adapt_mass_type;

  p | ppm_diffusion;
  p | ppm_flattening;
  p | ppm_minimum_pressure_support_parameter;
  p | ppm_pressure_free;
  p | ppm_steepening;
  p | ppm_use_minimum_pressure_support;

  p | field_uniform_density;

  p | physics_cosmology;
  p | physics_cosmology_hubble_constant_now;
  p | physics_cosmology_omega_lamda_now;
  p | physics_cosmology_omega_matter_now;
  p | physics_cosmology_omega_baryon_now;
  p | physics_cosmology_omega_cdm_now;
  p | physics_cosmology_comoving_box_size;
  p | physics_cosmology_max_expansion_rate;
  p | physics_cosmology_initial_redshift;
  p | physics_cosmology_final_redshift;

  p | physics_fluid_props_de_config;
  ::pup(p, physics_fluid_props_eos_variant);
  p | physics_fluid_props_fluid_floor_config;
  p | physics_fluid_props_mol_weight;

  p | physics_gravity;

  p | initial_bcenter_update_etot;

  p | initial_cloud_subsample_n;
  p | initial_cloud_radius;
  p | initial_cloud_center_x;
  p | initial_cloud_center_y;
  p | initial_cloud_center_z;
  p | initial_cloud_density_cloud;
  p | initial_cloud_density_wind;
  p | initial_cloud_velocity_wind;
  p | initial_cloud_etot_wind;
  p | initial_cloud_eint_wind;
  p | initial_cloud_metal_mass_frac;
  p | initial_cloud_initialize_uniform_bfield;
  PUParray(p,initial_cloud_uniform_bfield,3);
  p | initial_cloud_perturb_stddev;
  p | initial_cloud_trunc_dev;
  p | initial_cloud_perturb_seed;

  p | initial_cosmology_temperature;

  p | initial_collapse_rank;
  PUParray(p,initial_collapse_array,3);
  p | initial_collapse_radius_relative;
  p | initial_collapse_particle_ratio;
  p | initial_collapse_mass;
  p | initial_collapse_temperature;

  p | initial_grackle_test_minimum_H_number_density;
  p | initial_grackle_test_maximum_H_number_density;
  p | initial_grackle_test_minimum_temperature;
  p | initial_grackle_test_maximum_temperature;
  p | initial_grackle_test_minimum_metallicity;
  p | initial_grackle_test_maximum_metallicity;
  p | initial_grackle_test_reset_energies;

  p | initial_inclinedwave_alpha;
  p | initial_inclinedwave_beta;
  p | initial_inclinedwave_amplitude;
  p | initial_inclinedwave_lambda;
  p | initial_inclinedwave_parallel_vel;
  p | initial_inclinedwave_positive_vel;
  p | initial_inclinedwave_wave_type;

  p | initial_sedov_rank;
  PUParray(p,initial_sedov_array,3);
  p | initial_sedov_radius_relative;
  p | initial_sedov_pressure_in;
  p | initial_sedov_pressure_out;
  p | initial_sedov_density;

  PUParray(p,initial_sedov_random_array,3);
  p | initial_sedov_random_half_empty;
  p | initial_sedov_random_grackle_cooling;
  p | initial_sedov_random_max_blasts;
  p | initial_sedov_random_radius_relative;
  p | initial_sedov_random_pressure_in;
  p | initial_sedov_random_pressure_out;
  p | initial_sedov_random_density;
  p | initial_sedov_random_te_multiplier;

  p | initial_turbulence_density;
  p | initial_turbulence_pressure;
  p | initial_turbulence_temperature;

  p | initial_hdf5_max_level;
  p | initial_hdf5_format;
  PUParray(p, initial_hdf5_blocking,3);
  p | initial_hdf5_monitor_iter;
  p | initial_hdf5_field_files;
  p | initial_hdf5_field_datasets;
  p | initial_hdf5_field_names;
  p | initial_hdf5_field_coords;
  p | initial_hdf5_particle_files;
  p | initial_hdf5_particle_datasets;
  p | initial_hdf5_particle_coords;
  p | initial_hdf5_particle_types;
  p | initial_hdf5_particle_attributes;

  p | initial_music_field_coords;
  p | initial_music_field_datasets;
  p | initial_music_field_files;
  p | initial_music_field_names;
  p | initial_music_particle_attributes;
  p | initial_music_particle_coords;
  p | initial_music_particle_datasets;
  p | initial_music_particle_files;
  p | initial_music_particle_types;
  p | initial_music_throttle_close_count;
  p | initial_music_throttle_group_size;
  p | initial_music_throttle_internode;
  p | initial_music_throttle_intranode;
  p | initial_music_throttle_node_files;
  p | initial_music_throttle_seconds_delay;
  p | initial_music_throttle_seconds_stagger;

  p | initial_pm_field;
  p | initial_pm_mpp;
  p | initial_pm_level;

  p | initial_burkertbodenheimer_rank;
  PUParray(p,initial_burkertbodenheimer_array,3);
  p | initial_burkertbodenheimer_densityprofile;
  p | initial_burkertbodenheimer_mass;
  p | initial_burkertbodenheimer_outer_velocity;
  p | initial_burkertbodenheimer_particle_ratio;
  p | initial_burkertbodenheimer_radius_relative;
  p | initial_burkertbodenheimer_rotating;
  p | initial_burkertbodenheimer_temperature;

  PUParray(p, initial_feedback_test_position,3);
  p | initial_feedback_test_luminosity;
  p | initial_feedback_test_density;
  p | initial_feedback_test_e_density;
  p | initial_feedback_test_from_file;
  p | initial_feedback_test_HeI_density;
  p | initial_feedback_test_HeII_density;
  p | initial_feedback_test_HeIII_density;
  p | initial_feedback_test_HI_density;
  p | initial_feedback_test_HII_density;
  p | initial_feedback_test_metal_fraction;
  p | initial_feedback_test_star_mass;
  p | initial_feedback_test_temperature;

  PUParray(p, initial_IG_center_position,3);
  PUParray(p, initial_IG_bfield,3);
  p | initial_IG_analytic_velocity;
  p | initial_IG_disk_mass;
  p | initial_IG_disk_metal_fraction;
  p | initial_IG_disk_temperature;
  p | initial_IG_gas_fraction;
  p | initial_IG_gas_halo_density;
  p | initial_IG_gas_halo_mass;
  p | initial_IG_gas_halo_metal_fraction;
  p | initial_IG_gas_halo_radius;
  p | initial_IG_gas_halo_temperature;
  p | initial_IG_include_recent_SF;
  p | initial_IG_live_dm_halo;
  p | initial_IG_recent_SF_bin_size;
  p | initial_IG_recent_SF_end;
  p | initial_IG_recent_SF_seed;
  p | initial_IG_recent_SF_SFR;
  p | initial_IG_recent_SF_start;
  p | initial_IG_scale_height;
  p | initial_IG_scale_length;
  p | initial_IG_stellar_bulge;
  p | initial_IG_stellar_disk;
  p | initial_IG_use_gas_particles;

  p | initial_shock_tube_setup_name;
  p | initial_shock_tube_aligned_ax;
  p | initial_shock_tube_axis_velocity;
  p | initial_shock_tube_trans_velocity;
  p | initial_shock_tube_flip_initialize;

  p | initial_soup_rank;
  p | initial_soup_file;
  p | initial_soup_rotate;
  PUParray(p,initial_soup_array,3);
  PUParray(p,initial_soup_d_pos,3);
  PUParray(p,initial_soup_d_size,3);
  p | initial_soup_pressure_in;
  p | initial_soup_pressure_out;
  p | initial_soup_density;

  p | initial_merge_sinks_test_particle_data_filename;

  p | method_check_num_files;
  p | method_check_ordering;
  p | method_check_dir;
  p | method_check_monitor_iter;

  PUParray(p,initial_accretion_test_sink_position,3);
  PUParray(p,initial_accretion_test_sink_velocity,3);
  p | initial_accretion_test_sink_mass;
  p | initial_accretion_test_gas_density;
  p | initial_accretion_test_gas_pressure;
  p | initial_accretion_test_gas_radial_velocity;

  PUParray(p,initial_shu_collapse_center,3);
  PUParray(p,initial_shu_collapse_drift_velocity,3);
  p | initial_shu_collapse_truncation_radius;
  p | initial_shu_collapse_nominal_sound_speed;
  p | initial_shu_collapse_instability_parameter;
  p | initial_shu_collapse_external_density;
  p | initial_shu_collapse_central_sink_exists;
  p | initial_shu_collapse_central_sink_mass;

  PUParray(p,initial_bb_test_center,3);
  PUParray(p,initial_bb_test_drift_velocity,3);
  p | initial_bb_test_mean_density;
  p | initial_bb_test_fluctuation_amplitude;
  p | initial_bb_test_truncation_radius;
  p | initial_bb_test_nominal_sound_speed;
  p | initial_bb_test_angular_rotation_velocity;
  p | initial_bb_test_external_density;

  p | method_heat_alpha;

  p | method_hydro_method;
  p | method_hydro_dual_energy;
  p | method_hydro_dual_energy_eta_1;
  p | method_hydro_dual_energy_eta_2;
  p | method_hydro_reconstruct_method;
  p | method_hydro_reconstruct_conservative;
  p | method_hydro_reconstruct_positive;
  p | method_hydro_riemann_solver;

  p | method_feedback_flavor;
  p | method_feedback_ejecta_mass;
  p | method_feedback_supernova_energy;
  p | method_feedback_ejecta_metal_fraction;
  p | method_feedback_stencil;
  p | method_feedback_radius;
  p | method_feedback_shift_cell_center;
  p | method_feedback_ke_fraction;
  p | method_feedback_use_ionization_feedback;
  p | method_feedback_time_first_sn;

  p | method_feedback_supernovae;
  p | method_feedback_unrestricted_sn;
  p | method_feedback_stellar_winds;
  p | method_feedback_min_level;
  p | method_feedback_analytic_SNR_shell_mass;
  p | method_feedback_fade_SNR;
  p | method_feedback_NEvents;
  p | method_feedback_radiation;

  p | method_star_maker_flavor;
  p | method_star_maker_use_altAlpha;
  p | method_star_maker_use_density_threshold;
  p | method_star_maker_use_overdensity_threshold;
  p | method_star_maker_use_temperature_threshold;
  p | method_star_maker_use_critical_metallicity;
  p | method_star_maker_use_velocity_divergence;
  p | method_star_maker_use_dynamical_time;
  p | method_star_maker_use_cooling_time;
  p | method_star_maker_use_self_gravitating;
  p | method_star_maker_use_h2_self_shielding;
  p | method_star_maker_use_jeans_mass;
  p | method_star_maker_number_density_threshold;
  p | method_star_maker_overdensity_threshold;
  p | method_star_maker_critical_metallicity;
  p | method_star_maker_temperature_threshold;
  p | method_star_maker_maximum_mass_fraction;
  p | method_star_maker_efficiency;
  p | method_star_maker_minimum_star_mass;
  p | method_star_maker_maximum_star_mass;
  p | method_star_maker_min_level;
  p | method_star_maker_turn_off_probability;

  p | method_m1_closure;
  p | method_m1_closure_N_groups;
  p | method_m1_closure_flux_function;
  p | method_m1_closure_hll_file;
  p | method_m1_closure_clight_frac;
  p | method_m1_closure_photon_escape_fraction;
  p | method_m1_closure_radiation_spectrum;
  p | method_m1_closure_temperature_blackbody;
  p | method_m1_closure_particle_luminosity;
  p | method_m1_closure_SED;
  p | method_m1_closure_min_photon_density;
  p | method_m1_closure_attenuation;
  p | method_m1_closure_thermochemistry;
  p | method_m1_closure_recombination_radiation;
  p | method_m1_closure_H2_photodissociation;
  p | method_m1_closure_lyman_werner_background;
  p | method_m1_closure_LWB_J21;
  p | method_m1_closure_cross_section_calculator;
  p | method_m1_closure_sigmaN;
  p | method_m1_closure_sigmaE;
  p | method_m1_closure_energy_lower;
  p | method_m1_closure_energy_upper;
  p | method_m1_closure_energy_mean;

  p | method_turbulence_edot;

  p | method_gravity_grav_const;
  p | method_gravity_solver;
  p | method_gravity_order;
  p | method_gravity_dt_max;
  p | method_gravity_accumulate;

  p | method_background_acceleration_flavor;
  p | method_background_acceleration_mass;
  p | method_background_acceleration_DM_mass;
  p | method_background_acceleration_bulge_mass;
  p | method_background_acceleration_core_radius;
  p | method_background_acceleration_bulge_radius;
  p | method_background_acceleration_stellar_mass;
  p | method_background_acceleration_DM_mass_radius;
  p | method_background_acceleration_stellar_scale_height_r;
  p | method_background_acceleration_stellar_scale_height_z;
  p | method_background_acceleration_apply_acceleration;
  PUParray(p,method_background_acceleration_angular_momentum,3);
  PUParray(p,method_background_acceleration_center,3);

  p | method_pm_deposit_alpha;
  p | method_pm_update_max_dt;

  p | method_vlct_riemann_solver;
  p | method_vlct_half_dt_reconstruct_method;
  p | method_vlct_full_dt_reconstruct_method;
  p | method_vlct_theta_limiter;
  p | method_vlct_mhd_choice;

  p | method_merge_sinks_merging_radius_cells;

  p | method_accretion_accretion_radius_cells;
  p | method_accretion_flavor;
  p | method_accretion_physical_density_threshold_cgs;
  p | method_accretion_max_mass_fraction;

  p | method_sink_maker_jeans_length_resolution_cells,
  p | method_sink_maker_physical_density_threshold_cgs,
  p | method_sink_maker_check_density_maximum,
  p | method_sink_maker_max_mass_fraction,
  p | method_sink_maker_min_sink_mass_solar,
  p | method_sink_maker_max_offset_cell_fraction,
  p | method_sink_maker_offset_seed_shift,

  p | prolong_enzo_type;
  p | prolong_enzo_positive;
  p | prolong_enzo_use_linear;

  p | solver_pre_smooth;
  p | solver_post_smooth;
  p | solver_last_smooth;
  p | solver_coarse_solve;
  p | solver_domain_solve;
  p | solver_weight;
  p | solver_restart_cycle;
  p | solver_precondition;
  p | solver_coarse_level;
  p | solver_is_unigrid;

  p | stopping_redshift;

  p | units_mass;
  p | units_density;
  p | units_length;
  p | units_time;

  p | method_grackle_use_grackle;

  if (method_grackle_use_grackle) {
    p  | method_grackle_use_cooling_timestep;
    p  | method_grackle_radiation_redshift;
    p  | method_grackle_chemistry;
  }

}

//----------------------------------------------------------------------

void EnzoConfig::read(Parameters * p) throw()
{
  TRACE("BEGIN EnzoConfig::read()");

  // Read Cello parameters


  TRACE("EnzoCharm::read calling Config::read()");

  ((Config*)this) -> read (p);

  read_adapt_(p);

  read_field_(p);

  // Initial [sorted]
  read_initial_accretion_test_(p);
  read_initial_bb_test_(p);
  read_initial_bcenter_(p);
  read_initial_burkertbodenheimer_(p);
  read_initial_cloud_(p);
  read_initial_collapse_(p);
  read_initial_cosmology_(p);
  read_initial_feedback_test_(p);
  read_initial_grackle_(p);
  read_initial_hdf5_(p);
  read_initial_inclined_wave_(p);
  read_initial_isolated_galaxy_(p);
  read_initial_merge_sinks_test_(p);
  read_initial_music_(p);
  read_initial_pm_(p);
  read_initial_sedov_(p);
  read_initial_sedov_random_(p);
  read_initial_shock_tube_(p);
  read_initial_shu_collapse_(p);
  read_initial_soup_(p);
  read_initial_turbulence_(p);

  // it's important for read_physics_ to precede read_method_grackle_
  read_physics_(p);

  // Method [sorted]

  read_method_accretion_(p);
  read_method_background_acceleration_(p);
  read_method_check_(p);
  read_method_feedback_(p);
  read_method_grackle_(p);
  read_method_gravity_(p);
  read_method_heat_(p);
  read_method_merge_sinks_(p);
  read_method_pm_deposit_(p);
  read_method_pm_update_(p);
  read_method_ppm_(p);
  read_method_m1_closure_(p);
  read_method_sink_maker_(p);
  read_method_star_maker_(p);
  read_method_turbulence_(p);
  read_method_vlct_(p);

  read_prolong_enzo_(p);

  read_solvers_(p);

  read_stopping_(p);


  TRACE("END   EnzoConfig::read()");
}

//======================================================================

void EnzoConfig::read_adapt_(Parameters *p)
{

  adapt_mass_type.resize(num_adapt);

  for (int ia=0; ia<num_adapt; ia++) {

    std::string prefix = "Adapt:" + adapt_list[ia] + ":";
    adapt_mass_type[ia] = p->value_string(prefix+"mass_type","unknown");
    ASSERT2("EnzoConfig::read()",
	    "Unknown mass_type %s for parameter %s",
	    adapt_mass_type[ia].c_str(),(prefix+"mass_type").c_str(),
	    (adapt_type[ia] != "mass" ||
	     (adapt_mass_type[ia]=="dark" ||
	      adapt_mass_type[ia]=="baryon")));
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_field_(Parameters *p)
{
  field_uniform_density = p->value_float ("Field:uniform_density",1.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_collapse_(Parameters * p)
{
  initial_collapse_rank =  p->value_integer("Initial:collapse:rank",0);
  for (int i=0; i<initial_collapse_rank; i++) {
    initial_collapse_array[i] =
      p->list_value_integer (i,"Initial:collapse:array",1);
  }
  for (int i=initial_collapse_rank; i<3; i++) {
    initial_collapse_array[i] = 1;
  }
  initial_collapse_radius_relative =
    p->value_float("Initial:collapse:radius_relative",0.1);
  initial_collapse_particle_ratio =
    p->value_float("Initial:collapse:particle_ratio",0.0);
  initial_collapse_mass =
    p->value_float("Initial:collapse:mass",enzo_constants::mass_solar);
  initial_collapse_temperature =
    p->value_float("Initial:collapse:temperature",10.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_cosmology_(Parameters * p)
{
  initial_cosmology_temperature =
    p->value_float("Initial:cosmology:temperature",0.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_grackle_(Parameters * p)
{
  // Grackle test initialization
  initial_grackle_test_minimum_H_number_density =
    p->value_float("Initial:grackle_test:minimum_H_number_density",0.1);
  initial_grackle_test_maximum_H_number_density =
    p->value_float("Initial:grackle_test:maximum_H_number_density",1000.0);
  initial_grackle_test_minimum_temperature =
    p->value_float("Initial:grackle_test:minimum_temperature",10.0);
  initial_grackle_test_maximum_temperature =
    p->value_float("Initial:grackle_test:maximum_temperature",1.0E8);
  initial_grackle_test_minimum_metallicity =
    p->value_float("Initial:grackle_test:minimum_metallicity", 1.0E-4);
  initial_grackle_test_maximum_metallicity =
    p->value_float("Initial:grackle_test:maximum_metallicity", 1.0);
  initial_grackle_test_reset_energies =
    p->value_integer("Initial:grackle_test:reset_energies",0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_hdf5_(Parameters * p)
{
  const std::string name_initial = "Initial:hdf5:";

  initial_hdf5_max_level = p->value_integer (name_initial + "max_level", 0);
  initial_hdf5_format    = p->value_string  (name_initial + "format", "music");

  for (int i=0; i<3; i++) {
    initial_hdf5_blocking[i] =
      p->list_value_integer(i,name_initial+"blocking",1);
  }

  initial_hdf5_monitor_iter = p->value_integer (name_initial + "monitor_iter", 0);

  const int num_files = p->list_length (name_initial + "file_list");

  for (int index_file=0; index_file<num_files; index_file++) {

    std::string file_id = name_initial +
      p->list_value_string (index_file,name_initial+"file_list") + ":";

    const std::string type    = p->value_string (file_id + "type","");
    const std::string name    = p->value_string (file_id + "name","");
    const std::string file    = p->value_string (file_id + "file","");
    const std::string dataset = p->value_string (file_id + "dataset","");
    const std::string coords  = p->value_string (file_id + "coords","xyz");

    if (type == "particle") {

      const std::string attribute = p->value_string (file_id+"attribute","");

      initial_hdf5_particle_files.     push_back(file);
      initial_hdf5_particle_datasets.  push_back(dataset);
      initial_hdf5_particle_coords.    push_back(coords);
      initial_hdf5_particle_types.     push_back(name);
      initial_hdf5_particle_attributes.push_back(attribute);

    } else if (type == "field") {

      initial_hdf5_field_files.        push_back(file);
      initial_hdf5_field_datasets.     push_back(dataset);
      initial_hdf5_field_names.        push_back(name);
      initial_hdf5_field_coords.       push_back(coords);

    } else {
      ERROR2 ("EnzoConfig::read",
	      "Unknown particle type %s for parameter %s",
	      type.c_str(),(file_id+"type").c_str());
    }
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_music_(Parameters * p)
{
  const std::string name_initial = "Initial:music:";

  const int num_files = p->list_length (name_initial + "file_list");

  for (int index_file=0; index_file<num_files; index_file++) {

    std::string file_id = name_initial +
      p->list_value_string (index_file,name_initial+"file_list") + ":";

    std::string type    = p->value_string (file_id+"type","");
    std::string name    = p->value_string (file_id+"name","");
    std::string file    = p->value_string (file_id+"file","");
    std::string dataset = p->value_string (file_id+"dataset","");
    std::string coords  = p->value_string (file_id+"coords","xyz");

    if (type == "particle") {
      std::string attribute = p->value_string (file_id+"attribute","");
      //      if (name != "") {
      initial_music_particle_files.     push_back(file);
      initial_music_particle_datasets.  push_back(dataset);
      initial_music_particle_coords.    push_back(coords);
      initial_music_particle_types.     push_back(name);
      initial_music_particle_attributes.push_back(attribute);
      //      }
    } else if (type == "field") {

      initial_music_field_files.        push_back(file);
      initial_music_field_datasets.     push_back(dataset);
      initial_music_field_names.        push_back(name);
      initial_music_field_coords.       push_back(coords);
    } else {
      ERROR2 ("EnzoConfig::read",
	      "Unknown particle type %s for parameter %s",
	      type.c_str(),(file_id+"type").c_str());
    }
  }
  // "sleep_by_process", "limit_per_node"
  initial_music_throttle_internode = p->value_logical
    ("Initial:music:throttle_internode",false);
  initial_music_throttle_intranode = p->value_logical
    ("Initial:music:throttle_intranode",false);
  initial_music_throttle_node_files = p->value_logical
    ("Initial:music:throttle_node_files",false);
  initial_music_throttle_close_count = p->value_integer
    ("Initial:music:throttle_close_count",0);
  initial_music_throttle_group_size = p->value_integer
    ("Initial:music:throttle_group_size",std::numeric_limits<int>::max());
  initial_music_throttle_seconds_stagger = p->value_float
    ("Initial:music:throttle_seconds_stagger",0.0);
  initial_music_throttle_seconds_delay = p->value_float
    ("Initial:music:throttle_seconds_delay",0.0);

}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_pm_(Parameters * p)
{
  initial_pm_field        = p->value_string  ("Initial:pm:field","density");
  initial_pm_mpp          = p->value_float   ("Initial:pm:mpp",-1.0);
  initial_pm_level        = p->value_integer ("Initial:pm:level",-1);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_burkertbodenheimer_(Parameters * p)
{
  // Burkert Bodenheimer initialization

  initial_burkertbodenheimer_rank =  p->value_integer("Initial:burkertbodenheimer:rank",0);
  for (int i=0; i<initial_burkertbodenheimer_rank; i++) {
    initial_burkertbodenheimer_array[i] =
      p->list_value_integer (i,"Initial:burkertbodenheimer:array",1);
  }
  for (int i=initial_burkertbodenheimer_rank; i<3; i++) {
    initial_burkertbodenheimer_array[i] = 1;
  }
  initial_burkertbodenheimer_radius_relative =
    p->value_float("Initial:burkertbodenheimer:radius_relative",0.1);
  initial_burkertbodenheimer_particle_ratio =
    p->value_float("Initial:burkertbodenheimer:particle_ratio",0.0);
  initial_burkertbodenheimer_mass =
    p->value_float("Initial:burkertbodenheimer:mass",enzo_constants::mass_solar);
  initial_burkertbodenheimer_temperature =
    p->value_float("Initial:burkertbodenheimer:temperature",10.0);
  initial_burkertbodenheimer_densityprofile =
    p->value_integer ("Initial:burkertbodenheimer:densityprofile",2);
  initial_burkertbodenheimer_rotating =
   p->value_logical ("Initial:burkertbodenheimer:rotating",true);
  initial_burkertbodenheimer_outer_velocity =
   p->value_float ("Initial:burkertbodenheimer:outer_velocity",-1.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_inclined_wave_(Parameters * p)
{

  // InitialInclinedWave initialization

  initial_inclinedwave_alpha          = p->value_float
    ("Initial:inclined_wave:alpha",0.0);
  initial_inclinedwave_beta           = p->value_float
    ("Initial:inclined_wave:beta",0.0);
  initial_inclinedwave_amplitude      = p->value_float
    ("Initial:inclined_wave:amplitude",1.e-6);
  initial_inclinedwave_lambda         = p->value_float
    ("Initial:inclined_wave:lambda",1.0);
  // The default vaue for parallel_vel is known by EnzoInitialInclinedWave
  // to mean that a value was not specified
  initial_inclinedwave_parallel_vel   = p->value_float
    ("Initial:inclined_wave:parallel_vel", std::numeric_limits<double>::min());
  initial_inclinedwave_positive_vel   = p->value_logical
    ("Initial:inclined_wave:positive_vel",true);
  initial_inclinedwave_wave_type      = p->value_string
    ("Initial:inclined_wave:wave_type","alfven");
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_sedov_(Parameters * p)
{
  initial_sedov_rank = p->value_integer ("Initial:sedov:rank",0);

  initial_sedov_array[0] = p->list_value_integer (0,"Initial:sedov:array",1);
  initial_sedov_array[1] = p->list_value_integer (1,"Initial:sedov:array",1);
  initial_sedov_array[2] = p->list_value_integer (2,"Initial:sedov:array",1);

  initial_sedov_radius_relative =
    p->value_float("Initial:sedov:radius_relative",0.1);
  initial_sedov_pressure_in =
    p->value_float("Initial:sedov:pressure_in",1.0);
  initial_sedov_pressure_out =
    p->value_float("Initial:sedov:pressure_out",1e-5);
  initial_sedov_density =
    p->value_float("Initial:sedov:density",1.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_sedov_random_(Parameters * p)
{
  initial_sedov_random_array[0] =
    p->list_value_integer (0,"Initial:sedov_random:array",1);
  initial_sedov_random_array[1] =
    p->list_value_integer (1,"Initial:sedov_random:array",1);
  initial_sedov_random_array[2] =
    p->list_value_integer (2,"Initial:sedov_random:array",1);

  initial_sedov_random_half_empty =
    p->value_logical ("Initial:sedov_random:half_empty",false);
  initial_sedov_random_grackle_cooling =
    p->value_logical ("Initial:sedov_random:grackle_cooling",false);
  initial_sedov_random_max_blasts =
    p->value_integer ("Initial:sedov_random:max_blasts",1);
  initial_sedov_random_radius_relative =
    p->value_float   ("Initial:sedov_random:radius_relative",0.1);
  initial_sedov_random_pressure_in =
    p->value_float   ("Initial:sedov_random:pressure_in",1.0);
  initial_sedov_random_pressure_out =
    p->value_float   ("Initial:sedov_random:pressure_out",1e-5);
  initial_sedov_random_density =
    p->value_float   ("Initial:sedov_random:density",1.0);
  initial_sedov_random_te_multiplier =
    p->value_integer  ("Initial:sedov_random:te_multiplier",1);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_shock_tube_(Parameters * p)
{
  // Shock Tube Initialization
  initial_shock_tube_setup_name = p->value_string
    ("Initial:shock_tube:setup_name","");
  initial_shock_tube_aligned_ax = p->value_string
    ("Initial:shock_tube:aligned_ax","x");
  initial_shock_tube_axis_velocity = p->value_float
    ("Initial:shock_tube:axis_velocity",0.0);
  initial_shock_tube_trans_velocity = p->value_float
    ("Initial:shock_tube:transverse_velocity",0.0);
  initial_shock_tube_flip_initialize = p -> value_logical
    ("Initial:shock_tube:flip_initialize", false);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_bcenter_(Parameters * p)
{
  // VL+CT b-field initialization
  initial_bcenter_update_etot = p->value_logical
    ("Initial:vlct_bfield:update_etot",false);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_cloud_(Parameters * p)
{
  // Cloud Crush Initialization
  initial_cloud_subsample_n     = p->value_integer
    ("Initial:cloud:subsample_n",0);
  initial_cloud_radius          = p->value_float
    ("Initial:cloud:cloud_radius",0.0);
  initial_cloud_center_x        = p->value_float
    ("Initial:cloud:cloud_center_x",0.0);
  initial_cloud_center_y        = p->value_float
    ("Initial:cloud:cloud_center_y",0.0);
  initial_cloud_center_z        = p->value_float
    ("Initial:cloud:cloud_center_z",0.0);
  initial_cloud_density_cloud   = p->value_float
    ("Initial:cloud:cloud_density",0.0);
  initial_cloud_density_wind    = p->value_float
    ("Initial:cloud:wind_density",0.0);
  initial_cloud_velocity_wind   = p->value_float
    ("Initial:cloud:wind_velocity",0.0);
  initial_cloud_etot_wind       = p->value_float
    ("Initial:cloud:wind_total_energy",0.0);
  initial_cloud_eint_wind       = p->value_float
    ("Initial:cloud:wind_internal_energy",0.0);
  initial_cloud_metal_mass_frac = p->value_float
    ("Initial:cloud:metal_mass_fraction",0.0);
  initial_cloud_perturb_stddev  = p->value_float
    ("Initial:cloud:perturb_standard_deviation",0.0);
  initial_cloud_trunc_dev       = p->value_float
    ("Initial:cloud:perturb_truncation_deviation",0.0);
  int init_cloud_perturb_seed_  = p->value_integer
    ("Initial:cloud:perturb_seed",0);
  ASSERT("EnzoConfig::read()", "Initial:cloud:perturb_seed must be >=0",
	 init_cloud_perturb_seed_ >= 0);
  initial_cloud_perturb_seed = (unsigned int) init_cloud_perturb_seed_;

  int initial_cloud_uniform_bfield_length = p->list_length
    ("Initial:cloud:uniform_bfield");
  if (initial_cloud_uniform_bfield_length == 0){
    initial_cloud_initialize_uniform_bfield = false;
  } else if (initial_cloud_uniform_bfield_length == 3){
    initial_cloud_initialize_uniform_bfield = true;
    for (int i = 0; i <3; i++){
      initial_cloud_uniform_bfield[i] = p->list_value_float
	(i,"Initial:cloud:uniform_bfield");
    }
  } else {
    ERROR("EnzoConfig::read",
	  "Initial:cloud:uniform_bfield must contain 0 or 3 entries.");
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_soup_(Parameters * p)
{
  // InitialSoup initialization

  initial_soup_rank      = p->value_integer ("Initial:soup:rank",0);
  initial_soup_file      = p->value_string ("Initial:soup:file","soup.png");
  initial_soup_rotate    = p->value_logical ("Initial:soup:rotate",false);
  for (int axis=0; axis<3; axis++) {
    initial_soup_array[axis]  = p->list_value_integer
      (axis,"Initial:soup:array",1);
    initial_soup_d_pos[axis]  = p->list_value_float
      (axis,"Initial:soup:d_pos",0.0);
    initial_soup_d_size[axis] = p->list_value_float
      (axis,"Initial:soup:d_size",0.0);
  }
  initial_soup_pressure_in =
    p->value_float("Initial:soup:pressure_in",1.0);
  initial_soup_pressure_out =
    p->value_float("Initial:soup:pressure_out",1e-5);
  initial_soup_density =
    p->value_float("Initial:soup:density",1.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_turbulence_(Parameters * p)
{
  initial_turbulence_density = p->value_float
    ("Initial:turbulence:density",1.0);

  // Must specify pressure or temperature
  initial_turbulence_pressure =    p->value_float
    ("Initial:turbulence:pressure",   0.0);
  initial_turbulence_temperature = p->value_float
    ("Initial:turbulence:temperature",0.0);

  bool uses_turbulence = false;
  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "turbulence") uses_turbulence=true;
  }

  if (uses_turbulence) {
    ASSERT ("EnzoConfig::read",
  	    "Either initial turbulence pressure or temperature must be defined",
  	    ! ((initial_turbulence_pressure == 0.0) &&
  	       (initial_turbulence_temperature == 0.0)));
    ASSERT ("EnzoConfig::read",
  	    "Initial turbulence pressure and temperature cannot "
	    "both be defined",
  	    ! ((initial_turbulence_pressure != 0.0) &&
  	       (initial_turbulence_temperature != 0.0)));
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_isolated_galaxy_(Parameters * p)
{
  initial_IG_scale_length = p->value_float
    ("Initial:isolated_galaxy:scale_length", 0.0343218);
  initial_IG_scale_height = p->value_float
    ("Initial:isolated_galaxy:scale_height", 0.00343218);
  initial_IG_disk_mass = p->value_float
    ("Initial:isolated_galaxy:disk_mass", 42.9661);
  initial_IG_gas_fraction = p->value_float
    ("Initial:isolated_galaxy:gas_fraction", 0.2);
  initial_IG_disk_temperature = p->value_float
    ("Initial:isolated_galaxy:disk_temperature", 1.0E4);
  initial_IG_disk_metal_fraction = p->value_float
    ("Initial:isolated_galaxy:disk_metal_fraction", 1.0E-10);
  initial_IG_gas_halo_mass = p->value_float
    ("Initial:isolated_galaxy:gas_halo_mass", 0.1);
  initial_IG_gas_halo_temperature = p->value_float
    ("Initial:isolated_galaxy:gas_halo_temperature", 1.0E4);
  initial_IG_gas_halo_density = p->value_float
    ("Initial:isolated_galaxy:gas_halo_density", 0.0);
  initial_IG_gas_halo_radius = p->value_float
    ("Initial:isolated_galaxy:gas_halo_radius", 1.0);
  initial_IG_gas_halo_metal_fraction = p->value_float
    ("Initial:isolated_galaxy:gas_halo_metal_fraction", 1.0E-10);
  initial_IG_use_gas_particles = p->value_logical
    ("Initial:isolated_galaxy:use_gas_particles", false);
  initial_IG_live_dm_halo = p->value_logical
    ("Initial:isolated_galaxy:live_dm_halo",false);
  initial_IG_stellar_disk = p->value_logical
    ("Initial:isolated_galaxy:stellar_disk", false);
  initial_IG_stellar_bulge = p->value_logical
    ("Initial:isolated_galaxy:stellar_bulge", false);
  initial_IG_analytic_velocity = p->value_logical
    ("Initial:isolated_galaxy:analytic_velocity", false);
  initial_IG_include_recent_SF = p->value_logical
    ("Initial:isolated_galaxy:include_recent_SF", false);
  initial_IG_recent_SF_start = p->value_float
    ("Initial:isolated_galaxy:recent_SF_start", -100.0);
  initial_IG_recent_SF_end = p->value_float
    ("Initial:isolated_galaxy:recent_SF_end", 0.0);
  initial_IG_recent_SF_SFR = p->value_float
    ("Initial:isolated_galaxy:recent_SF_SFR", 2.0);
  initial_IG_recent_SF_bin_size = p->value_float
    ("Initial:isolated_galaxy:recent_SF_bin_size", 5.0);
  initial_IG_recent_SF_seed = p->value_integer
    ("Initial:isolated_galaxy:recent_SF_seed", 12345);

  for (int axis=0; axis<3; axis++) {
    initial_IG_center_position[axis]  = p->list_value_float
      (axis,"Initial:isolated_galaxy:center_position",0.5);
    initial_IG_bfield[axis] = p->list_value_float
      (axis, "Initial:isolated_galaxy:bfield",0.0);
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_initial_feedback_test_(Parameters * p)
{
  for (int axis=0; axis<3; axis++){
    initial_feedback_test_position[axis] = p->list_value_float
      (axis, "Initial:feedback_test:position", 0.5);
  }
  initial_feedback_test_luminosity = p->value_float
    ("Initial:feedback_test:luminosity", 0.0);

  initial_feedback_test_density = p->value_float
    ("Initial:feedback_test:density", 1.0E-24);

  initial_feedback_test_HI_density = p->value_float
    ("Initial:feedback_test:HI_density", 1.0E-24);

  initial_feedback_test_HII_density = p->value_float
    ("Initial:feedback_test:HII_density", 1.0E-100);

  initial_feedback_test_HeI_density = p->value_float
    ("Initial:feedback_test:HeI_density", 1.0E-100);

  initial_feedback_test_HeII_density = p->value_float
    ("Initial:feedback_test:HeII_density", 1.0E-100);

  initial_feedback_test_HeIII_density = p->value_float
    ("Initial:feedback_test:HeIII_density", 1.0E-100);

  initial_feedback_test_e_density = p->value_float
    ("Initial:feedback_test:e_density", 1.0E-100);

  initial_feedback_test_star_mass = p->value_float
    ("Initial:feedback_test:star_mass", 1000.0);

  initial_feedback_test_temperature = p->value_float
    ("Initial:feedback_test:temperature", 1.0E4);

  initial_feedback_test_from_file = p->value_logical
    ("Initial:feedback_test:from_file", false);

  initial_feedback_test_metal_fraction = p->value_float
    ("Initial:feedback_test:metal_fraction", 0.01);
}

void EnzoConfig::read_initial_merge_sinks_test_(Parameters * p)
{
  initial_merge_sinks_test_particle_data_filename= p->value_string
    ("Initial:merge_sinks_test:particle_data_filename","");
}

void EnzoConfig::read_initial_accretion_test_(Parameters * p)
{
  for (int axis=0; axis<3; axis++){
    initial_accretion_test_sink_position[axis] = p->list_value_float
      (axis, "Initial:accretion_test:sink_position", 0.0);
  }

  for (int axis=0; axis<3; axis++){
    initial_accretion_test_sink_velocity[axis] = p->list_value_float
      (axis, "Initial:accretion_test:sink_velocity", 0.0);
  }

  initial_accretion_test_sink_mass = p->value_float
    ("Initial:accretion_test:sink_mass",0.0);

  initial_accretion_test_gas_density = p->value_float
    ("Initial:accretion_test:gas_density",1.0e-6);

  initial_accretion_test_gas_pressure = p->value_float
    ("Initial:accretion_test:gas_pressure",1.0e-6);

  initial_accretion_test_gas_radial_velocity = p->value_float
    ("Initial:accretion_test:gas_radial_velocity",0.0);
}

void EnzoConfig::read_initial_shu_collapse_(Parameters * p)
{
  for (int axis=0; axis<3; axis++){
    initial_shu_collapse_center[axis] = p->list_value_float
      (axis, "Initial:shu_collapse:center", 0.0);
  }

  for (int axis=0; axis<3; axis++){
    initial_shu_collapse_drift_velocity[axis] = p->list_value_float
      (axis, "Initial:shu_collapse:drift_velocity", 0.0);
  }

  initial_shu_collapse_truncation_radius = p->value_float
    ("Initial:shu_collapse:truncation_radius",1.0);

  initial_shu_collapse_nominal_sound_speed = p->value_float
    ("Initial:shu_collapse:nominal_sound_speed",1.0);

  initial_shu_collapse_instability_parameter = p->value_float
    ("Initial:shu_collapse:instability_parameter",2.1);

  initial_shu_collapse_external_density = p->value_float
    ("Initial:shu_collapse:external_density",1.0e-6);

  initial_shu_collapse_central_sink_exists = p->value_logical
    ("Initial:shu_collapse:central_sink_exists",false);

  initial_shu_collapse_central_sink_mass = p->value_float
    ("Initial:shu_collapse:central_sink_mass",0.0);
}

void EnzoConfig::read_initial_bb_test_(Parameters * p)
{
  for (int axis=0; axis<3; axis++){
    initial_bb_test_center[axis] = p->list_value_float
      (axis, "Initial:bb_test:center", 0.0);
  }

  for (int axis=0; axis<3; axis++){
    initial_bb_test_drift_velocity[axis] = p->list_value_float
      (axis, "Initial:bb_test:drift_velocity", 0.0);
  }

  initial_bb_test_mean_density = p->value_float
    ("Initial:bb_test:mean_density",1.0e-6);

  initial_bb_test_fluctuation_amplitude = p->value_float
    ("Initial:bb_test:fluctuation_amplitude",0.0);

  initial_bb_test_truncation_radius = p->value_float
    ("Initial:bb_test:truncation_radius",1.0);

  initial_bb_test_nominal_sound_speed = p->value_float
    ("Initial:bb_test:nominal_sound_speed",1.0);

  initial_bb_test_angular_rotation_velocity = p->value_float
    ("Initial:bb_test:angular_rotation_velocity",0.0);

  initial_bb_test_external_density = p->value_float
    ("Initial:bb_test:external_density",1.0e-6);
}

void EnzoConfig::read_method_grackle_(Parameters * p)

{
  method_grackle_use_grackle = false;

  /// Grackle parameters

  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "grackle") method_grackle_use_grackle=true;
  }

  // Defaults alert PUP::er() to ignore
  if (method_grackle_use_grackle) {

    method_grackle_use_cooling_timestep = p->value_logical
      ("Method:grackle:use_cooling_timestep", false);

    // for when not using cosmology - redshift of UVB
    method_grackle_radiation_redshift = p->value_float
      ("Method:grackle:radiation_redshift", -1.0);

    // Now, we will initialize method_grackle_chemistry
    // - we do this with a factory method that directly examines the parameter
    //   values within the "Method:grackle:*" group.
    // - Because Grackle has so many parameter values, it's very easy to make a
    //   small mistake when specifying the name of a parameter value and not
    //   notice until much later. For that reason, the factory method
    //   aggressively reports unexpected parameters as errors.
    // - to help this method, we provide 2 sets of parameter names

    //   1. specify all of the Grackle parameters that we will manually setup
    //      based on the values passed for other Enzo-E parameters. Errors will
    //      be reported if any of these are encountered
    const std::unordered_set<std::string> forbid_leaf_names = {"use_grackle",
                                                               "Gamma"};

    //   2. specify all parameters that MAY occur within the "Method:grackle:*"
    //      group that should be ignored by the factory method. (This needs to
    //      be updated if we introduce additional parameters for configuring
    //      EnzoMethodGrackle)
    const std::unordered_set<std::string> ignore_leaf_names =
      {"use_cooling_timestep", "radiation_redshift",
       // the next option is deprecated and is only listed in the short-term
       // for backwards compatability (it should now be replaced by
       // "Physics:fluid_props:floors:metallicity")
       "metallicity_floor",
       // for backwards compatability, we manually use "data_file" to
       // initialize "grackle_data_file" parameter (in the future, we may want
       // to change this)
       "data_file", "grackle_data_file",
       // the final two parameters auto-parsed by other Cello machinery
       "type", "courant"};

    method_grackle_chemistry = GrackleChemistryData::from_parameters
      (*p, "Method:grackle", forbid_leaf_names, ignore_leaf_names);

    // now let's manually initialize the handful of remaining runtime
    // parameters that are stored within method_grackle_chemistry

    // 1. use "Method:grackle:data_file" to initialize "grackle_data_file" 
    if (p->param("Method:grackle:grackle_data_file") != nullptr){
      ERROR("EnzoConfig::read_method_grackle_",
            "for backwards compatability, the user can't specify "
            "\"Method:grackle:grackle_data_file\". Instead, they must specify "
            "\"Method:grackle:data_file\".");
    } else if (p->param("Method:grackle:data_file") != nullptr) {
      std::string fname = p->value_string("Method:grackle:data_file", "");
      ASSERT("EnzoConfig::read_method_grackle_",
             "\"Method:grackle:data_file\" can't be an empty string",
             fname.length() > 0); // sanity check!
      method_grackle_chemistry.set<std::string>("grackle_data_file", fname);
    } else {
      ERROR("EnzoConfig::read_method_grackle_",
            "\"Method:grackle:data_file\" is required when using grackle");
    }

    // 2. update the value of use_grackle
    method_grackle_chemistry.set<int>("use_grackle",
                                      method_grackle_use_grackle);

    // 3. Copy over parameters from Enzo-E to Grackle
    if (physics_fluid_props_eos_variant.holds_alternative<EnzoEOSIdeal>()) {
      method_grackle_chemistry.set<double>
        ("Gamma", physics_fluid_props_eos_variant.get<EnzoEOSIdeal>().gamma);
    } else {
      ERROR("EnzoConfig::read_method_grackle_",
            "Grackle currently can't be used when Enzo-E is configured to use "
            "an equation of state other than the ideal gas law");
    }

    // In the future, we may want to manually set use_radiative_transfer based
    // on an Enzo-E parameter for turning RT on / off:
    //method_grackle_chemistry.set<int>("use_radiative_transfer", ENZO_E_PARAMETER_NAME);
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_feedback_(Parameters * p)
{
  method_feedback_flavor = p->value_string
    ("Method:feedback:flavor","distributed");

  method_feedback_ejecta_mass = p->value_float
    ("Method:feedback:ejecta_mass",0.0);

  method_feedback_supernova_energy = p->value_float
    ("Method:feedback:supernova_energy",1.0);

  method_feedback_ejecta_metal_fraction = p->value_float
    ("Method:feedback:ejecta_metal_fraction",0.1);

  method_feedback_stencil = p->value_integer
    ("Method:feedback:stencil",3);

  method_feedback_radius = p->value_float
    ("Method:feedback:radius",-1.0);

  method_feedback_shift_cell_center = p->value_logical
    ("Method:feedback:shift_cell_center", true);

  method_feedback_ke_fraction = p->value_float
    ("Method:feedback:ke_fraction", 0.0);

  method_feedback_time_first_sn = p->value_float
    ("Method:feedback:time_first_sn", -1.0);

  method_feedback_use_ionization_feedback = p->value_logical
    ("Method:feedback:use_ionization_feedback", false);

  // MethodFeedbackSTARSS parameters
  method_feedback_supernovae = p->value_logical
    ("Method:feedback:supernovae",true);

  method_feedback_unrestricted_sn = p->value_logical
    ("Method:feedback:unrestricted_sn",true);

  method_feedback_stellar_winds = p->value_logical
    ("Method:feedback:stellar_winds",true);

  method_feedback_min_level = p->value_integer
    ("Method:feedback:min_level",0);

  method_feedback_analytic_SNR_shell_mass = p->value_logical
    ("Method:feedback:analytic_SNR_shell_mass",true);

  method_feedback_fade_SNR = p->value_logical
    ("Method:feedback:fade_SNR",true);

  method_feedback_NEvents = p->value_integer
    ("Method:feedback:NEvents",-1);

  method_feedback_radiation = p->value_logical
    ("Method:feedback:radiation", true);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_star_maker_(Parameters * p)
{
  method_star_maker_flavor = p->value_string
    ("Method:star_maker:flavor","stochastic");

  method_star_maker_use_altAlpha = p->value_logical
    ("Method:star_maker:use_altAlpha",false);

  method_star_maker_use_density_threshold = p->value_logical
    ("Method:star_maker:use_density_threshold",false);

  method_star_maker_use_overdensity_threshold = p->value_logical
    ("Method:star_maker:use_overdensity_threshold",false);

  method_star_maker_use_velocity_divergence = p->value_logical
    ("Method:star_maker:use_velocity_divergence",false);

  method_star_maker_use_dynamical_time = p->value_logical
    ("Method:star_maker:use_dynamical_time",false);

  method_star_maker_use_cooling_time = p->value_logical
    ("Method:star_maker:use_cooling_time",false);

  method_star_maker_use_self_gravitating = p->value_logical
    ("Method:star_maker:use_self_gravitating", false);

  method_star_maker_use_h2_self_shielding = p->value_logical
    ("Method:star_maker:use_h2_self_shielding", false);

  method_star_maker_use_jeans_mass = p->value_logical
    ("Method:star_maker:use_jeans_mass", false);

  method_star_maker_use_temperature_threshold = p->value_logical
    ("Method:star_maker:use_temperature_threshold",false);

  method_star_maker_use_critical_metallicity = p->value_logical
    ("Method:star_maker:use_critical_metallicity",false);

  method_star_maker_number_density_threshold = p->value_float
    ("Method:star_maker:number_density_threshold",0.0);

  method_star_maker_overdensity_threshold = p->value_float
    ("Method:star_maker:overdensity_threshold",0.0);

  method_star_maker_temperature_threshold = p->value_float
    ("Method:star_maker:temperature_threshold",1.0E4);

  method_star_maker_critical_metallicity = p->value_float
    ("Method:star_maker:critical_metallicity",0.0);

  method_star_maker_maximum_mass_fraction = p->value_float
    ("Method:star_maker:maximum_mass_fraction",0.05);

  method_star_maker_efficiency = p->value_float
    ("Method:star_maker:efficiency",0.01);

  method_star_maker_minimum_star_mass = p->value_float
    ("Method:star_maker:minimum_star_mass",0.0);

  method_star_maker_maximum_star_mass = p->value_float
    ("Method:star_maker:maximum_star_mass",-1.0);

  method_star_maker_min_level = p->value_integer
    ("Method:star_maker:min_level",0);

  method_star_maker_turn_off_probability = p->value_logical
    ("Method:star_maker:turn_off_probability",false);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_m1_closure_(Parameters * p)
{
  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "m1_closure") method_m1_closure=true;
  }

  method_m1_closure_N_groups = p->value_integer
    ("Method:m1_closure:N_groups",1);

  method_m1_closure_flux_function = p->value_string
    ("Method:m1_closure:flux_function","GLF");

  method_m1_closure_hll_file = p->value_string
    ("Method:m1_closure:hll_file","hll_evals.list");

  method_m1_closure_clight_frac = p->value_float
    ("Method:m1_closure:clight_frac",1.0);

  method_m1_closure_photon_escape_fraction = p->value_float
    ("Method:m1_closure:photon_escape_fraction",1.0);

  method_m1_closure_radiation_spectrum = p->value_string
    ("Method:m1_closure:radiation_spectrum","custom");

  method_m1_closure_temperature_blackbody = p->value_float
    ("Method:m1_closure:temperature_blackbody",0.0);

  method_m1_closure_particle_luminosity = p->value_float
    ("Method:m1_closure:particle_luminosity",-1.0);

  method_m1_closure_min_photon_density = p->value_float
    ("Method:m1_closure:min_photon_density",0.0);

  method_m1_closure_attenuation = p->value_logical
    ("Method:m1_closure:attenuation", true);

  method_m1_closure_thermochemistry = p->value_logical
    ("Method:m1_closure:thermochemistry", true);

  method_m1_closure_recombination_radiation = p->value_logical
    ("Method:m1_closure:recombination_radiation",false);

  method_m1_closure_H2_photodissociation = p->value_logical
    ("Method:m1_closure:H2_photodissociation", false);

  method_m1_closure_lyman_werner_background = p->value_logical
    ("Method:m1_closure:lyman_werner_background", false);

  method_m1_closure_LWB_J21 = p->value_float
    ("Method:m1_closure:LWB_J21", -1.0);

  method_m1_closure_cross_section_calculator = p->value_string
    ("Method:m1_closure:cross_section_calculator","vernier");

  method_m1_closure_SED.resize(method_m1_closure_N_groups);
  method_m1_closure_energy_lower.resize(method_m1_closure_N_groups);
  method_m1_closure_energy_upper.resize(method_m1_closure_N_groups);
  method_m1_closure_energy_mean.resize(method_m1_closure_N_groups);
  
  int N_species = 3; // number of ionizable species for RT (HI, HeI, HeII)
  method_m1_closure_sigmaN.resize(method_m1_closure_N_groups * N_species);
  method_m1_closure_sigmaE.resize(method_m1_closure_N_groups * N_species);

  // make default energy bins equally spaced between 1 eV and 101 eV
  double bin_width = 100.0 / method_m1_closure_N_groups;
  for (int i=0; i<method_m1_closure_N_groups; i++) {
    // default SED (if this is being used) is flat spectrum
    method_m1_closure_SED[i] = p->list_value_float
      (i,"Method:m1_closure:SED", 1.0/method_m1_closure_N_groups);

    method_m1_closure_energy_lower[i] = p->list_value_float
      (i,"Method:m1_closure:energy_lower", 1.0 + bin_width*i);

    method_m1_closure_energy_upper[i] = p->list_value_float
      (i,"Method:m1_closure:energy_upper", 1.0 + bin_width*(i+1));

    method_m1_closure_energy_mean[i] = p->list_value_float
      (i,"Method:m1_closure:energy_mean", 
         0.5*(method_m1_closure_energy_lower[i] + method_m1_closure_energy_upper[i]));

    for (int j=0; j<N_species; j++) {
      int sig_index = i*N_species + j;
      method_m1_closure_sigmaN[sig_index] = p->list_value_float
        (sig_index,"Method:m1_closure:sigmaN", 0.0);

      method_m1_closure_sigmaE[sig_index] = p->list_value_float
        (sig_index,"Method:m1_closure:sigmaE", 0.0);
    }
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_background_acceleration_(Parameters * p)
{
  method_background_acceleration_flavor = p->value_string
   ("Method:background_acceleration:flavor","unknown");

  method_background_acceleration_mass = p->value_float
   ("Method:background_acceleration:mass",0.0);

  method_background_acceleration_DM_mass = p->value_float
   ("Method:background_acceleration:DM_mass",-1.0);

  method_background_acceleration_bulge_mass = p->value_float
    ("Method:background_acceleration:bulge_mass", 0.0);

  method_background_acceleration_core_radius = p->value_float
    ("Method:background_acceleration:core_radius", 1.0E-10);

  method_background_acceleration_bulge_radius = p->value_float
    ("Method:background_acceleration:bulge_radius", 1.0E-10);

  method_background_acceleration_stellar_mass = p->value_float
    ("Method:background_acceleration:stellar_mass", 0.0);

  method_background_acceleration_DM_mass_radius = p->value_float
   ("Method:background_acceleration:DM_mass_radius", 0.0);

  method_background_acceleration_stellar_scale_height_r = p->value_float
   ("Method:background_acceleration:stellar_scale_height_r", 1.0E-10);

  method_background_acceleration_stellar_scale_height_z = p->value_float
   ("Method:background_acceleration:stellar_scale_height_z", 1.0E-10);

  method_background_acceleration_apply_acceleration = p->value_logical
    ("Method:background_acceleration:apply_acceleration", true);

  for (int axis = 0; axis < 3; axis++){
    method_background_acceleration_center[axis] = p->list_value_float
      (axis,"Method:background_acceleration:center",0.5);
    method_background_acceleration_angular_momentum[axis] = p->list_value_float
      (axis,"Method:background_acceleration:angular_momentum",0);
  }

  // Not sure if I need. Seems this flag tells the hydo solver
  // if gravity exists... so I would expect to need this... but Does
  // not get triggered for self-gravity at the moment... so not sure
  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "background_acceleration") physics_gravity=true;
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_vlct_(Parameters * p)
{
  method_vlct_riemann_solver = p->value_string
    ("Method:mhd_vlct:riemann_solver","hlld");
  method_vlct_half_dt_reconstruct_method = p->value_string
    ("Method:mhd_vlct:half_dt_reconstruct_method","nn");
  method_vlct_full_dt_reconstruct_method = p->value_string
    ("Method:mhd_vlct:full_dt_reconstruct_method","plm");
  method_vlct_theta_limiter = p->value_float
    ("Method:mhd_vlct:theta_limiter", 1.5);

  // we should raise an error if mhd_choice is not specified
  bool uses_vlct = false;
  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "mhd_vlct") uses_vlct=true;
  }
  method_vlct_mhd_choice = p->value_string
    ("Method:mhd_vlct:mhd_choice", "");
  if (uses_vlct && (method_vlct_mhd_choice == "")){
    ERROR("EnzoConfig::read", "Method:mhd_vlct:mhd_choice was not specified");
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_gravity_(Parameters * p)
{
  method_gravity_grav_const = p->value_float
    ("Method:gravity:grav_const",6.67384e-8);

  method_gravity_solver = p->value_string
    ("Method:gravity:solver","unknown");

  //--------------------------------------------------
  // Physics
  //--------------------------------------------------

  method_gravity_order = p->value_integer
    ("Method:gravity:order",4);

  method_gravity_accumulate = p->value_logical
    ("Method:gravity:accumulate",true);

  method_gravity_dt_max = p->value_float
    ("Method:gravity:dt_max",1.0e10);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_check_(Parameters * p)
{
  p->group_set(0,"Method");
  p->group_push("check");

  method_check_num_files = p->value_integer
    ("num_files",1);
  method_check_ordering = p->value_string
    ("ordering","order_morton");

  if (p->type("dir") == parameter_string) {
    method_check_dir.resize(1);
    method_check_dir[0] = p->value_string("dir","");
  } else if (p->type("dir") == parameter_list) {
    int size = p->list_length("dir");
    if (size > 0) method_check_dir.resize(size);
    for (int i=0; i<size; i++) {
      method_check_dir[i] = p->list_value_string(i,"dir","");
    }
  }
  method_check_monitor_iter = p->value_integer("monitor_iter",0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_heat_(Parameters * p)
{
  method_heat_alpha = p->value_float
    ("Method:heat:alpha",1.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_merge_sinks_(Parameters * p)
{
  method_merge_sinks_merging_radius_cells = p->value_float
    ("Method:merge_sinks:merging_radius_cells",8.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_accretion_(Parameters * p)
{
  method_accretion_accretion_radius_cells = p->value_float
    ("Method:accretion:accretion_radius_cells",4.0);
  method_accretion_flavor = p->value_string
    ("Method:accretion:flavor","");
  method_accretion_physical_density_threshold_cgs = p->value_float
    ("Method:accretion:physical_density_threshold_cgs",1.0e-24);
  method_accretion_max_mass_fraction = p->value_float
    ("Method:accretion:max_mass_fraction",0.25);

}

//----------------------------------------------------------------------

void EnzoConfig::read_method_sink_maker_(Parameters * p)
{
  method_sink_maker_jeans_length_resolution_cells = p->value_float
    ("Method:sink_maker:jeans_length_resolution_cells",4.0);
  method_sink_maker_physical_density_threshold_cgs = p->value_float
    ("Method:sink_maker:physical_density_threshold_cgs",1.0e-24);
  method_sink_maker_check_density_maximum = p->value_logical
    ("Method:sink_maker:check_density_maximum",true);
  method_sink_maker_max_mass_fraction = p->value_float
    ("Method:sink_maker:max_mass_fraction",0.25);
  method_sink_maker_min_sink_mass_solar = p->value_float
    ("Method:sink_maker:min_sink_mass_solar",0.0);
  method_sink_maker_max_offset_cell_fraction = p->value_float
    ("Method:sink_maker:max_offset_cell_fraction",0.0);
  int method_sink_maker_offset_seed_shift_input = p->value_integer
    ("Method:sink_maker:offset_seed_shift",0);
  ASSERT("EnzoConfig::read()", "Method:sink_maker:offset_seed_shift must be >=0",
	 method_sink_maker_offset_seed_shift_input >= 0);
  method_sink_maker_offset_seed_shift = (uint64_t) method_sink_maker_offset_seed_shift_input;
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_pm_deposit_(Parameters * p)
{
  method_pm_deposit_alpha = p->value_float ("Method:pm_deposit:alpha",0.5);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_pm_update_(Parameters * p)
{
  method_pm_update_max_dt = p->value_float
    ("Method:pm_update:max_dt", std::numeric_limits<double>::max());
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_ppm_(Parameters * p)
{
  double floor_default = 1e-6;

  ppm_diffusion = p->value_logical
    ("Method:ppm:diffusion", false);
  ppm_flattening = p->value_integer
    ("Method:ppm:flattening", 3);
  ppm_minimum_pressure_support_parameter = p->value_integer
    ("Method:ppm:minimum_pressure_support_parameter",100);
  ppm_pressure_free = p->value_logical
    ("Method:ppm:pressure_free",false);
  ppm_steepening = p->value_logical
    ("Method:ppm:steepening", false);
  ppm_use_minimum_pressure_support = p->value_logical
    ("Method:ppm:use_minimum_pressure_support",false);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_turbulence_(Parameters * p)
{
  method_turbulence_edot = p->value_float
    ("Method:turbulence:edot",-1.0);
  method_turbulence_mach_number = p->value_float
    ("Method:turbulence:mach_number",0.0);
}

//----------------------------------------------------------------------

void EnzoConfig::read_physics_(Parameters * p)
{
  num_physics = p->list_length("Physics:list");

  for (int index_physics=0; index_physics<num_physics; index_physics++) {

    std::string name =
      p->list_value_string(index_physics,"Physics:list");

    std::string full_name = std::string("Physics:") + name;

    if (physics_list[index_physics] == "cosmology") {

      physics_cosmology = true;

      physics_cosmology_hubble_constant_now = p->value_float
        (full_name + ":hubble_constant_now",0.701);

      physics_cosmology_omega_matter_now = p->value_float
        (full_name + ":omega_matter_now",   0.279);

      physics_cosmology_omega_baryon_now = p->value_float
        (full_name + ":omega_baryon_now",   1.0);

      physics_cosmology_omega_cdm_now = p->value_float
        (full_name + ":omega_cdm_now",   0.0);

      physics_cosmology_omega_lamda_now = p->value_float
        (full_name + ":omega_lambda_now",   0.721);


      physics_cosmology_comoving_box_size = p->value_float
        (full_name + ":comoving_box_size", 64.0);

      physics_cosmology_max_expansion_rate = p->value_float
        (full_name + ":max_expansion_rate", 0.01);

      physics_cosmology_initial_redshift = p->value_float
        (full_name + ":initial_redshift",  20.0);;

      physics_cosmology_final_redshift = p->value_float
        (full_name + ":final_redshift",  0.0);;

    }

    if (physics_list[index_physics] == "fluid_prop"){
      ERROR("EnzoConfig::read_physics_",
            "\"fluid_prop\" is a typo for \"fluid_props\"");
    }

    if (physics_list[index_physics] == "gravity") {

      physics_gravity = true;

    }
  }

  // this is intentionally done outside of the for-loop (for
  // backwards-compatability purposes)
  read_physics_fluid_props_(p);
}

//----------------------------------------------------------------------

namespace{

  /// parse a parameter that is allowed to be a float or a list of floats
  ///
  /// returns an empty vector if the parameter does not exist
  std::vector<double> coerce_param_list_(Parameters * p,
                                         const std::string& parameter)
  {
    std::vector<double> out;
    if (p->type(parameter) == parameter_float) {
      out.push_back(p->value_float(parameter));
    } else if (p->type(parameter) == parameter_list) {
      const int list_length = p->list_length(parameter);
      for (int i = 0; i < list_length; i++){
        out.push_back(p->list_value_float(i, parameter));
      }
    } else if (p->param(parameter) != nullptr) {
      ERROR1("coerce_param_list_",
             "The \"%s\" parameter was specified with an invalid type. When "
             "specified, it must be a float or list of floats",
             parameter.c_str());
    }
    return out;
  }

  //----------------------------------------------------------------------

  EnzoDualEnergyConfig parse_de_config_(Parameters * p,
                                        const std::string& hydro_type)
  {
    EnzoDualEnergyConfig out = EnzoDualEnergyConfig::build_disabled();

    // fetch names of parameters in Physics:fluid_props:dual_energy
    p->group_set(0, "Physics");
    p->group_set(1, "fluid_props");
    p->group_set(2, "dual_energy");
    std::vector<std::string> names = p->leaf_parameter_names();

    const bool missing_de_config = names.size() == 0;
    if (!missing_de_config){ // parse Physics:fluid_props:dual_energy
      const std::string type = p->value_string
        ("Physics:fluid_props:dual_energy:type", "disabled");

      const std::string eta_paramname = "Physics:fluid_props:dual_energy:eta";
      const bool eta_exists = p->param(eta_paramname) != nullptr;
      const std::vector<double> eta_list = coerce_param_list_(p, eta_paramname);

      // raise an error if parameters were specified if there are unexpected
      // parameters. We are being a little extra careful here.
      for (const std::string& name : names){
        ASSERT1("parse_de_config_",
                "Unexpected parameter: \"Physics:fluid_props:dual_energy:%s\"",
                name.c_str(), (name == "type") | (name == "eta"));
      }

      // now actually construct the output object
      if (type == "disabled"){
        ASSERT1("parse_de_config_",
                "when dual energy is disabled, \"%s\" can't be specified",
                eta_paramname.c_str(), !eta_exists);
        out = EnzoDualEnergyConfig::build_disabled();
      } else if (type == "modern"){
        ASSERT3("parse_de_config_",
                "\"%s\" was used to specify %d values. When specified for the "
                "\"%s\" dual energy formalism, it must provide 1 value.",
                eta_paramname.c_str(), (int)names.size(), type.c_str(),
                (eta_list.size() == 1) | !eta_exists);
        double eta = eta_exists ? eta_list[0] : 0.001;
        out = EnzoDualEnergyConfig::build_modern_formulation(eta);
      } else if (type == "bryan95"){
        ASSERT3("parse_de_config_",
                "\"%s\" was used to specify %d value(s). When specified for "
                "the \"%s\" dual energy formalism, it must provide 2 value.",
                eta_paramname.c_str(), (int)names.size(), type.c_str(),
                (eta_list.size() == 2) | !eta_exists);
        double eta_1 = eta_exists ? eta_list[0] : 0.001;
        double eta_2 = eta_exists ? eta_list[1] : 0.1;
        out = EnzoDualEnergyConfig::build_bryan95_formulation(eta_1, eta_2);
      } else {
        ERROR1("parse_de_config_",
               "\"Physics:fluid_props:dual_energy:type\" is invalid: \"%s\"",
               type.c_str());
      }
    }

    // look for dual energy parameters specified within the hydro solver (for
    // backwards compatibility)
    if ((hydro_type != "") && (hydro_type != "ppml")) {
      std::string legacy_de_param = "Method:" + hydro_type + ":dual_energy";
      bool legacy_param_exists = p->param(legacy_de_param) != nullptr;

      if (legacy_param_exists & !missing_de_config){
        ERROR1("parse_de_config_",
               "legacy parameter, \"%s\", duplicates other parameters",
               legacy_de_param.c_str());
      } else if (legacy_param_exists) {
        WARNING1("parse_de_config_",
                 "\"%s\" is a legacy parameter that will be removed",
                 legacy_de_param.c_str());
        bool use_de = p->value_logical(legacy_de_param, false);
        if (!use_de){
          out = EnzoDualEnergyConfig::build_disabled();
        } else if (hydro_type == "ppm"){
          out = EnzoDualEnergyConfig::build_bryan95_formulation
            (p->value_float("Method:ppm:dual_energy_eta_1", 0.001),
             p->value_float("Method:ppm:dual_energy_eta_2", 0.1));
        } else {
          out = EnzoDualEnergyConfig::build_modern_formulation
            (p->value_float("Method:mhd_vlct:dual_energy_eta", 0.001));
        }
      }
    }
    return out;
  }

  //----------------------------------------------------------------------

  EnzoEOSVariant parse_eos_choice_(Parameters * p,
                                   const std::string& hydro_type)
  {
    // Prior to the creation of the EnzoEOSVariant class, Enzo-E effectively
    // assumed at a global level that an ideal EOS was in use and stored gamma
    // at a global level.
    //
    // - gamma's value was originally parsed from "Field:gamma" and later from
    //   "Physics:fluid_props:eos:gamma". As an aside, while it was always the
    //   plan to have the ``EnzoPhysicsFluidProps`` class track the EOS type,
    //   the ability to track EOS type was added a while after the fact (in the
    //   same PR that introduced ``EnzoEOSVariant``).
    // - gamma's default value has always been 5/3.
    // - when the Ppml solver was used, it simply ignored the value of gamma
    //   and internally used an Isothermal EOS.
    // - technically, extension points were put into place within the vl+ct
    //   solver to support other types of solvers, but those were never used.


    // get the name of EnzoEOSIdeal (useful since it's the default eos type)
    const std::string ideal_name = EnzoEOSIdeal::name();

    // check whether the legacy parameter was specified
    const bool legacy_gamma_specified = p->param("Field:gamma") != nullptr;

    // fetch names of parameters in Physics:fluid_props:eos
    p->group_set(0, "Physics");
    p->group_set(1, "fluid_props");
    p->group_set(2, "eos");
    std::vector<std::string> names = p->leaf_parameter_names();

    const bool missing_eos_config = names.size() == 0;

    if (legacy_gamma_specified && !missing_eos_config) {
      ERROR("parse_eos_choice_",
            "\"Field:gamma\" isn't valid since parameters are specified "
            "within the \"Physics:fluid_props:eos\" parameter group");

    } else if (!missing_eos_config) {
      // this branch does the main work of the function

      // STEP 1: define some useful variables
      const std::string prefix = "Physics:fluid_props:eos:";
      const bool is_type_specified = p->param(prefix + "type") != nullptr;
      // following variable is used for maintaining backwards compatability
      const bool is_gamma_specified = p->param(prefix + "gamma") != nullptr;

      // STEP 2: determine the eos-type
      std::string type;
      if (is_type_specified) {
        type = p->value(prefix + "type","");
      } else if (is_gamma_specified) {
        WARNING1("parse_eos_choice_",
                 "Going forward, \"Physics:fluid_props:eos:type\" must be set "
                 "when there are other parameters in the subgroup. For "
                 "backwards compatability, this is being set to \"%s\" since "
                 "the only other parameter in that group is \"gamma\"",
                 ideal_name.c_str());
        type = ideal_name;
      } else {
        ERROR("parse_eos_choice_",
              "\"Physics:fluid_props:eos:type\" must be set when there are "
              "other parameters in the subgroup.");
      }

      // STEP 3: actually build the EOS object and return it
      if (type == ideal_name) { // EnzoEOSIdeal

        // this case is a little funky, since we allow type to not actually be
        // a parameter (for backwards compatability purposes).
        std::size_t num_params = (1 + (std::size_t)(is_type_specified));
        ASSERT1("parse_eos_choice_",
                "the only allowed parameters are \"type\" and \"gamma\" in "
                "the \"Physics:fluid_props:eos\" parameter group when making "
                "an \"%s\" eos", type.c_str(),
                (num_params == names.size()) && is_gamma_specified);
        double gamma = p->value_float(prefix+"gamma", -1.0);
        return EnzoEOSVariant(EnzoEOSIdeal::construct(gamma));

      } else if ( type == EnzoEOSIsothermal::name() ){

        ASSERT1("parse_eos_choice_",
                "when building an \"%s\" eos, \"type\" is the only parameter "
                "allowed in the \"Physics:fluid_props:eos\" parameter group ",
                type.c_str(), (names.size() == 1) && is_type_specified);
        return EnzoEOSVariant(EnzoEOSIsothermal());

      } else {
        ERROR1("parse_eos_choice_",
               "there is currently no support for building of type \"%s\".",
               type.c_str());
      }
    }


    if (legacy_gamma_specified) {
      WARNING1("parse_eos_choice_",
               "\"Field:gamma\" is a legacy parameter that will be removed. "
               "It is being used to configure an \"%s\" EOS. Going forward, "
               "set parameters in the \"Physics:fluid_props:eos\" parameter "
               "group instead.",
               ideal_name.c_str());
      double gamma = p->value_float("Field:gamma", -1.0);
      return EnzoEOSVariant(EnzoEOSIdeal::construct(gamma));

    } else if (hydro_type == "ppml") {
      std::string type = EnzoEOSIsothermal::name();
      WARNING1("parse_eos_choice_",
               "Defaulting to \"%s\" EOS since for backwards compatability "
               "since the PPML solver is in use and no parameters were set "
               "in the \"Physics:fluid_props:eos\" parameter group. In the "
               "future, this behavior will be dropped.",
               type.c_str());
      return EnzoEOSVariant(EnzoEOSIsothermal());

    } else {
      const double default_gamma = 5.0/3.0;
      WARNING2("parse_eos_choice_",
               "No parameters specified in the \"Physics:fluid_props:eos\" "
               "parameter group. Defaulting to an \"%s\" eos with gamma = "
               "%#.16g.",
               ideal_name.c_str(), default_gamma);
      return EnzoEOSVariant(EnzoEOSIdeal::construct(default_gamma));

    }

  }

  //----------------------------------------------------------------------

  EnzoFluidFloorConfig parse_fluid_floor_config_(Parameters * p,
                                                 const std::string& hydro_type,
                                                 bool using_grackle)
  {
    // initialize default values (a value <= 0 means there is no floor)
    double density_floor = 0.0;
    double pressure_floor = 0.0;
    double temperature_floor = 0.0;
    double metal_mass_frac_floor = 0.0;

    auto get_ptr_to_floor_var = [&](const std::string name)
      {
        if (name == "density") { return &density_floor; }
        if (name == "pressure") { return &pressure_floor; }
        if (name == "temperature") { return &temperature_floor; }
        if (name == "metallicity") { return &metal_mass_frac_floor; }
        return (double*)nullptr;
      };

    // fetch names of parameters in Physics:fluid_props:floors. If any of them
    // exist, let's parse them
    p->group_set(0, "Physics");
    p->group_set(1, "fluid_props");
    p->group_set(2, "floors");
    std::vector<std::string> floor_l = p->leaf_parameter_names();

    const bool no_legacy = (floor_l.size() > 0);
    if (no_legacy){
      for (const std::string& name : floor_l){
        double* ptr = get_ptr_to_floor_var(name);
        if (ptr == nullptr){
          ERROR1("EnzoConfig::read_physics_fluid_props_",
                 "no support for placing a floor on \"%s\"", name.c_str());
        } else if (name == "metallicity") {
          *ptr = (p->value_float(p->full_name(name)) *
                  enzo_constants::metallicity_solar);
        } else {
          *ptr = p->value_float(p->full_name(name));
        }
      }
    }

    // now let's consider the legacy options (for the appropriate hydro solver
    // and Grackle). if there were no parameters in Physics:fluid_props:floors,
    // let's parse them. Otherwise, let's raise an error
    const std::vector<std::array<std::string,3>> legacy_params =
      {{"density", "ppm", "density_floor"},
       {"pressure", "ppm", "pressure_floor"},
       {"temperature", "ppm", "temperature_floor"},
       {"density", "mhd_vlct", "density_floor"},
       {"pressure", "mhd_vlct", "pressure_floor"},
       {"metallicity", "grackle", "metallicity_floor"}};

    for (const std::array<std::string,3>& triple : legacy_params){
      if ( (("grackle" != triple[1]) && (hydro_type != triple[1])) ||
           (("grackle" == triple[1]) && (!using_grackle)) ) {
        continue;
      }
      std::string full_name = "Method:" + triple[1] + ":" + triple[2];
      if (p->param(full_name) == nullptr) {continue;}
      if (no_legacy){
        ERROR1("EnzoConfig::read_physics_fluid_props_",
               "legacy parameter \"%s\" is invalid since the "
               "\"Physics:fluid_props:floors\" parameters are specified",
               full_name.c_str());
      } else {
        WARNING2("EnzoConfig::read_physics_fluid_props_",
                 "\"%s\" is a deprecated parameter (it will be removed in the "
                 "future). Use \"Physics:fluid_props:floors:%s\" instead.",
                 full_name.c_str(), triple[0].c_str());
        if (triple[0] == "metallicity"){
          *(get_ptr_to_floor_var(triple[0]))
            = p->value_float(full_name) * enzo_constants::metallicity_solar;
        } else {
          *(get_ptr_to_floor_var(triple[0])) = p->value_float(full_name);
        }
      }
    }

    return {density_floor, pressure_floor, temperature_floor,
            metal_mass_frac_floor};
  }

}

//----------------------------------------------------------------------

void EnzoConfig::read_physics_fluid_props_(Parameters * p)
{
  // determine the hydro method (if any) so we know which legacy parameters to
  // look for.
  const std::vector<std::string>& mlist = this->method_list;
  bool has_ppm = std::find(mlist.begin(), mlist.end(), "ppm") != mlist.end();
  bool has_ppml = std::find(mlist.begin(), mlist.end(), "ppml") != mlist.end();
  bool has_vlct = std::find(mlist.begin(), mlist.end(),
                            "mhd_vlct") != mlist.end();
  std::string hydro_type = "";
  if ((int(has_ppm) + int(has_ppml) + int(has_vlct)) > 1){
    ERROR("EnzoConfig::read_physics_fluid_props_",
          "a given simulation can only use up to 1 of the following solvers: "
          "{\"ppm\", \"ppml\", \"mhd_vlct\"}");
  } else if (has_ppm){
    hydro_type = "ppm";
  } else if (has_ppml){
    hydro_type = "ppml";
  } else if (has_vlct){
    hydro_type = "mhd_vlct";
  }
  bool has_grackle = std::find(mlist.begin(), mlist.end(),
                               "grackle") != mlist.end();

  // determine the dual energy formalism configuration
  physics_fluid_props_de_config = parse_de_config_(p, hydro_type);

  // determine the fluid floor configuration
  physics_fluid_props_fluid_floor_config =
    parse_fluid_floor_config_(p, hydro_type, has_grackle);

  // determine the nominal choice of the EOS (the EOS is currently independent
  // of the molecular weight)
  physics_fluid_props_eos_variant = parse_eos_choice_(p, hydro_type);

  // determine molecular weight
  {
    double default_val = 0.6;
    double legacy_value = p->value_float("Method:ppm:mol_weight", -1);
    double actual_value = p->value_float("Physics:fluid_props:mol_weight", -1);

    if (legacy_value == -1) {
      if (actual_value == -1) { actual_value = default_val; }
      physics_fluid_props_mol_weight = actual_value;
    } else if (actual_value == -1) {
      WARNING("EnzoConfig::read_physics_fluid_props_",
              "\"Method:ppm:mol_weight\" is a legacy parameter that will be "
              "removed.");
      physics_fluid_props_mol_weight = legacy_value;
    } else {
      ERROR("EnzoConfig::read_physics_fluid_props_",
            "\"Method:ppm:mol_weight\" isn't valid since "
            "\"Physics:fluid_props:mol_weight\" is specified.");
    }
  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_prolong_enzo_(Parameters * p)
{
  prolong_enzo_type       = p->value_string  ("Prolong:enzo:type","2A");
  prolong_enzo_positive   = p->value_logical ("Prolong:enzo:positive",true);
  prolong_enzo_use_linear = p->value_logical ("Prolong:enzo:use_linear",false);
}

//----------------------------------------------------------------------

void EnzoConfig::read_solvers_(Parameters * p)
{
  num_solvers = p->list_length("Solver:list");

  solver_pre_smooth.  resize(num_solvers);
  solver_coarse_solve.resize(num_solvers);
  solver_domain_solve.resize(num_solvers);
  solver_post_smooth. resize(num_solvers);
  solver_last_smooth. resize(num_solvers);
  solver_weight.      resize(num_solvers);
  solver_restart_cycle.resize(num_solvers);
  solver_precondition.resize(num_solvers);
  solver_coarse_level.resize(num_solvers);
  solver_is_unigrid.resize(num_solvers);

  for (int index_solver=0; index_solver<num_solvers; index_solver++) {

    std::string solver_name =
      std::string("Solver:") + p->list_value_string(index_solver,"Solver:list");

    std::string solver;

    solver = p->value_string (solver_name + ":precondition","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_precondition[index_solver] = solver_index[solver];
    } else {
      solver_precondition[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":pre_smooth","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_pre_smooth[index_solver] = solver_index[solver];
    } else {
      solver_pre_smooth[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":coarse_solve","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_coarse_solve[index_solver] = solver_index[solver];
    } else {
      solver_coarse_solve[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":domain_solve","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_domain_solve[index_solver] = solver_index[solver];
    } else {
      solver_domain_solve[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":post_smooth","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_post_smooth[index_solver] = solver_index[solver];
    } else {
      solver_post_smooth[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":last_smooth","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_last_smooth[index_solver] = solver_index[solver];
    } else {
      solver_last_smooth[index_solver] = -1;
    }

    solver_weight[index_solver] =
      p->value_float(solver_name + ":weight",1.0);

    solver_restart_cycle[index_solver] =
      p->value_integer(solver_name + ":restart_cycle",1);

    solver_coarse_level[index_solver] =
      p->value_integer (solver_name + ":coarse_level",
                        solver_min_level[index_solver]);

    solver_is_unigrid[index_solver] =
      p->value_logical (solver_name + ":is_unigrid",false);

  }
}

//----------------------------------------------------------------------

void EnzoConfig::read_stopping_(Parameters * p)
{
  stopping_redshift = p->value_float ("Stopping:redshift",0.0);
}

//----------------------------------------------------------------------
