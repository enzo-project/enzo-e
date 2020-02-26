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

#ifdef CONFIG_USE_GRACKLE
// Operator to allow Grackle's chemistry data to PUP
inline void operator|(PUP::er &p, chemistry_data &c){
 // all values are single ints, floats, or doubles with the
 // exception of grackle_data_file
 p | c.use_grackle;
 p | c.with_radiative_cooling;
 p | c.primordial_chemistry;
 p | c.metal_cooling;
 p | c.UVbackground;

 int length = (c.grackle_data_file == NULL) ? 0 : strlen(c.grackle_data_file);
 p | length;
 if (length > 0){
   if (p.isUnpacking()){
     c.grackle_data_file=new char[length+1];
   }
   PUParray(p, c.grackle_data_file,length+1);
 } else {
   c.grackle_data_file = NULL;
 }

 p | c.cmb_temperature_floor;
 p | c.Gamma;
 p | c.h2_on_dust;
 p | c.photoelectric_heating;
 p | c.photoelectric_heating_rate;
 p | c.use_volumetric_heating_rate;
 p | c.use_specific_heating_rate;
 p | c.three_body_rate;
 p | c.cie_cooling;
 p | c.h2_optical_depth_approximation;
 p | c.ih2co;
 p | c.ipiht;
 p | c.HydrogenFractionByMass;
 p | c.DeuteriumToHydrogenRatio;
 p | c.SolarMetalFractionByMass;
 p | c.NumberOfTemperatureBins;
 p | c.CaseBRecombination;
 p | c.TemperatureStart;
 p | c.TemperatureEnd;
 p | c.NumberOfDustTemperatureBins;
 p | c.DustTemperatureStart;
 p | c.DustTemperatureEnd;
 p | c.Compton_xray_heating;
 p | c.LWbackground_sawtooth_suppression;
 p | c.LWbackground_intensity;
 p | c.UVbackground_redshift_on;
 p | c.UVbackground_redshift_off;
 p | c.UVbackground_redshift_fullon;
 p | c.UVbackground_redshift_drop;
 p | c.cloudy_electron_fraction_factor;
 p | c.use_radiative_transfer;
 p | c.radiative_transfer_coupled_rate_solver;
 p | c.radiative_transfer_intermediate_step;
 p | c.radiative_transfer_hydrogen_only;
 p | c.self_shielding_method;
 p | c.H2_self_shielding;
}
#endif

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
      physics_gravity(false),
      // EnzoInitialCosmology
      initial_cosmology_temperature(0.0),
      // EnzoInitialCollapse
      initial_collapse_rank(0),
      initial_collapse_radius_relative(0.0),
      initial_collapse_particle_ratio(0.0),
      initial_collapse_mass(0.0),
      initial_collapse_temperature(0.0),
      // EnzoGrackleTest
#ifdef CONFIG_USE_GRACKLE
      initial_grackle_test_maximum_H_number_density(1000.0),
      initial_grackle_test_maximum_metallicity(1.0),
      initial_grackle_test_maximum_temperature(1.0E8),
      initial_grackle_test_minimum_H_number_density(0.1),
      initial_grackle_test_minimum_metallicity(1.0E-4),
      initial_grackle_test_minimum_temperature(10.0),
      initial_grackle_test_reset_energies(0),
#endif /* CONFIG_USE_GRACKLE */
      initial_feedback_test_density(),
      initial_feedback_test_star_mass(),
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
      // EnzoInitialIsolatedGalaxy
      //   AE: Maybe these values (and those in cpp) don't matter
      //       are they overwritten by the read-in (even when not found in param file)?
      initial_IG_scale_length(0.0343218),       // Gas disk scale length in code units
      initial_IG_scale_height(0.00343218),      // Gas disk scale height in code units
      initial_IG_disk_mass(42.9661),            // Gas disk mass in code units
      initial_IG_gas_fraction(0.2),             // Gas disk M_gas / M_star
      initial_IG_disk_temperature(1e4),         // Gas disk temperature in K
      initial_IG_disk_metal_fraction(1.0E-10),         // Gas disk metal fraction
      initial_IG_gas_halo_mass(0.1),             // Gas halo total mass in code units
      initial_IG_gas_halo_temperature(1e4),      // Gas halo initial temperature
      initial_IG_gas_halo_metal_fraction(1.0E-10),      // Gas halo metal fraction
      initial_IG_gas_halo_density(0.0),          // Gas halo uniform density (ignored if zero)
      initial_IG_gas_halo_radius(1.0),           // Gas halo maximum radius in code units
      initial_IG_use_gas_particles(false),       //
      initial_IG_live_dm_halo(false),
      initial_IG_stellar_disk(false),
      initial_IG_stellar_bulge(false),
      initial_IG_analytic_velocity(false),
      initial_IG_include_recent_SF(false),
      initial_IG_recent_SF_start(-100.0),
      initial_IG_recent_SF_end(0.0),
      initial_IG_recent_SF_bin_size(5.0),
      initial_IG_recent_SF_SFR(2.0),
      initial_IG_recent_SF_seed(12345),
      // EnzoProlong
      interpolation_method(""),
      // EnzoMethodCheckGravity
      method_check_gravity_particle_type(),
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
      /// EnzoMethodFeedback
      method_feedback_ejecta_mass(0.0),
      method_feedback_supernova_energy(1.0),
      method_feedback_ejecta_metal_fraction(0.0),
      method_feedback_stencil(3),
      method_feedback_shift_cell_center(true),
      method_feedback_ke_fraction(0.0),
      method_feedback_use_ionization_feedback(false),
      /// EnzoMethodStarMaker
      method_star_maker_type(""),
      method_star_maker_use_density_threshold(true),           // check above density threshold before SF
      method_star_maker_use_velocity_divergence(true),         // check for converging flow before SF
      method_star_maker_use_dynamical_time(true),              //
      method_star_maker_use_self_gravitating(false),           // 
      method_star_maker_use_h2_self_shielding(false),
      method_star_maker_use_jeans_mass(false),
      method_star_maker_number_density_threshold(0.0),      // Number density threshold in cgs
      method_star_maker_maximum_mass_fraction(0.5),            // maximum cell mass fraction to convert to stars
      method_star_maker_efficiency(0.01),            // star maker efficiency
      method_star_maker_minimum_star_mass(1.0E4),    // minium star particle mass in solar masses
      method_star_maker_maximum_star_mass(1.0E4),    // maximum star particle mass in solar masses
      // EnzoMethodNull
      method_null_dt(0.0),
      // EnzoMethodTurbulence
      method_turbulence_edot(0.0),
      method_turbulence_mach_number(0.0),
      // EnzoMethodGrackle
      method_grackle_use_grackle(false),
#ifdef CONFIG_USE_GRACKLE
      method_grackle_chemistry(),
      method_grackle_use_cooling_timestep(false),
      method_grackle_radiation_redshift(-1.0),
#endif
      // EnzoMethodGravity
      method_gravity_grav_const(0.0),
      method_gravity_solver(""),
      method_gravity_order(4),
      method_gravity_accumulate(false),
      // EnzoMethodBackgroundAcceleration
      method_background_acceleration_type(""),
      method_background_acceleration_mass(0.0),
      method_background_acceleration_DM_mass(0.0),
      method_background_acceleration_DM_density(0.0),
      method_background_acceleration_bulge_mass(0.0),
      method_background_acceleration_core_radius(0.0),
      method_background_acceleration_bulge_radius(0.0),
      method_background_acceleration_stellar_mass(0.0),
      method_background_acceleration_DM_mass_radius(0.0),
      method_background_acceleration_stellar_scale_height_r(0.0),
      method_background_acceleration_stellar_scale_height_z(0.0),
      method_background_acceleration_apply_acceleration(true),
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
      initial_IG_center_position[axis] = 0.5;
      initial_IG_bfield[axis]         = 0.0;
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

  /// EnzoGrackleTest
#ifdef CONFIG_USE_GRACKLE
  double                     initial_grackle_test_maximum_H_number_density;
  double                     initial_grackle_test_maximum_metallicity;
  double                     initial_grackle_test_maximum_temperature;
  double                     initial_grackle_test_minimum_H_number_density;
  double                     initial_grackle_test_minimum_metallicity;
  double                     initial_grackle_test_minimum_temperature;
  int                        initial_grackle_test_reset_energies;
#endif /* CONFIG_USE_GRACKLE */

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
  double                     initial_feedback_test_density;
  double                     initial_feedback_test_star_mass;

  /// EnzoInitialIsolatedGalaxy
  double                     initial_IG_center_position[3];
  double                     initial_IG_bfield[3];
  double                     initial_IG_scale_length;
  double                     initial_IG_scale_height;
  double                     initial_IG_disk_mass;
  double                     initial_IG_gas_fraction;
  double                     initial_IG_disk_temperature;
  double                     initial_IG_disk_metal_fraction;
  double                     initial_IG_gas_halo_mass;
  double                     initial_IG_gas_halo_temperature;
  double                     initial_IG_gas_halo_metal_fraction;
  double                     initial_IG_gas_halo_density;
  double                     initial_IG_gas_halo_radius;
  bool                       initial_IG_use_gas_particles;
  bool                       initial_IG_live_dm_halo;
  bool                       initial_IG_stellar_bulge;
  bool                       initial_IG_stellar_disk;
  bool                       initial_IG_analytic_velocity;
  bool                       initial_IG_include_recent_SF;
  double                     initial_IG_recent_SF_start;
  double                     initial_IG_recent_SF_end;
  double                     initial_IG_recent_SF_bin_size;
  double                     initial_IG_recent_SF_SFR;
  int                        initial_IG_recent_SF_seed;

  /// EnzoProlong
  std::string                interpolation_method;

  /// EnzoMethodHeat
  double                     method_heat_alpha;

  /// EnzoMethodCheckGravity
  std::string                method_check_gravity_particle_type;
  
  /// EnzoMethodHydro
  std::string                method_hydro_method;
  bool                       method_hydro_dual_energy;
  double                     method_hydro_dual_energy_eta_1;
  double                     method_hydro_dual_energy_eta_2;
  std::string                method_hydro_reconstruct_method;
  bool                       method_hydro_reconstruct_conservative;
  bool                       method_hydro_reconstruct_positive;
  std::string                method_hydro_riemann_solver;

  /// EnzoMethodFeedback

  double                    method_feedback_ejecta_mass;
  double                    method_feedback_supernova_energy;
  double                    method_feedback_ejecta_metal_fraction;
  double                    method_feedback_ke_fraction;
  int                       method_feedback_stencil;
  bool                      method_feedback_shift_cell_center;
  bool                      method_feedback_use_ionization_feedback;

  /// EnzoMethodStarMaker

  std::string               method_star_maker_type;
  bool                      method_star_maker_use_density_threshold;
  bool                      method_star_maker_use_velocity_divergence;
  bool                      method_star_maker_use_dynamical_time;
  bool                      method_star_maker_use_h2_self_shielding;
  bool                      method_star_maker_use_jeans_mass;
  bool                      method_star_maker_use_self_gravitating;
  double                    method_star_maker_number_density_threshold;
  double                    method_star_maker_maximum_mass_fraction;
  double                    method_star_maker_efficiency;
  double                    method_star_maker_minimum_star_mass;
  double                    method_star_maker_maximum_star_mass;



  /// EnzoMethodNull
  double                     method_null_dt;

  /// EnzoMethodTurbulence
  double                     method_turbulence_edot;
  double                     method_turbulence_mach_number;

  /// EnzoMethodGrackle
  bool                       method_grackle_use_grackle;
#ifdef CONFIG_USE_GRACKLE
  chemistry_data *           method_grackle_chemistry;
  bool                       method_grackle_use_cooling_timestep;
  double                     method_grackle_radiation_redshift;
#endif /* CONFIG_USE_GRACKLE */

  /// EnzoMethodGravity
  double                     method_gravity_grav_const;
  std::string                method_gravity_solver;
  int                        method_gravity_order;
  bool                       method_gravity_accumulate;

  /// EnzoMethodBackgroundAcceleration

  std::string                method_background_acceleration_type;
  double                     method_background_acceleration_mass;
  double                     method_background_acceleration_DM_mass;
  double                     method_background_acceleration_DM_density;
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

  /// Mg0 coarse grid solver

  std::vector<int>           solver_coarse_level;
  std::vector<int>           solver_is_unigrid;

  /// Stop at specified redshift for cosmology
  double                     stopping_redshift;

};

extern EnzoConfig g_enzo_config;

#endif /* PARAMETERS_ENZO_CONFIG_HPP */
