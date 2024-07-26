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
  // EnzoInitialCosmology
  initial_cosmology_temperature(0.0),
  // EnzoInitialCollapse
  initial_collapse_rank(0),
  initial_collapse_radius_relative(0.0),
  initial_collapse_particle_ratio(0.0),
  initial_collapse_mass(0.0),
  initial_collapse_temperature(0.0),
  // EnzoInitialHdf5
  initial_hdf5_max_level(),
  initial_hdf5_format(),
  initial_hdf5_blocking(),
  initial_hdf5_monitor_iter(),
  initial_hdf5_field_files(),
  initial_hdf5_field_datasets(),
  initial_hdf5_field_names(),
  initial_hdf5_field_coords(),
  initial_hdf5_field_levels(),
  initial_hdf5_particle_files(),
  initial_hdf5_particle_datasets(),
  initial_hdf5_particle_coords(),
  initial_hdf5_particle_types(),
  initial_hdf5_particle_attributes(),
  initial_hdf5_particle_levels(),
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
  method_check_include_ghosts(false),
  // EnzoInitialMergeSinksTest
  initial_merge_sinks_test_particle_data_filename(""),
  // EnzoInitialAccretionTest
  initial_accretion_test_sink_mass(0.0),
  initial_accretion_test_gas_density(0.0),
  initial_accretion_test_gas_pressure(0.0),
  initial_accretion_test_gas_radial_velocity(0.0),
  // EnzoInitialBBTest
  initial_bb_test_mean_density(0.0),
  initial_bb_test_fluctuation_amplitude(0.0),
  initial_bb_test_truncation_radius(0.0),
  initial_bb_test_nominal_sound_speed(0.0),
  initial_bb_test_angular_rotation_velocity(0.0),
  initial_bb_test_external_density(0.0),
  // EnzoMethodGravity
  method_gravity_type_super(),
  // EnzoMethodInference
  method_inference_level_base(0),
  method_inference_level_array(0),
  method_inference_level_infer(0),
  method_inference_field_group(),
  method_inference_overdensity_threshold(0),
  // EnzoMethodTurbulence
  method_turbulence_edot(0.0),
  method_turbulence_mach_number(0.0),
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
    initial_sedov_array[i] = 0;
    initial_collapse_array[i] = 0;
    initial_IG_center_position[i] = 0.5;
    initial_IG_bfield[i] = 0.0;

  }
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

  p | initial_bcenter_update_etot;

  p | initial_cosmology_temperature;

  p | initial_collapse_rank;
  PUParray(p,initial_collapse_array,3);
  p | initial_collapse_radius_relative;
  p | initial_collapse_particle_ratio;
  p | initial_collapse_mass;
  p | initial_collapse_temperature;

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
  p | initial_hdf5_field_levels;
  p | initial_hdf5_particle_files;
  p | initial_hdf5_particle_datasets;
  p | initial_hdf5_particle_coords;
  p | initial_hdf5_particle_types;
  p | initial_hdf5_particle_attributes;
  p | initial_hdf5_particle_levels;

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

  p | initial_merge_sinks_test_particle_data_filename;

  p | method_check_num_files;
  p | method_check_ordering;
  p | method_check_dir;
  p | method_check_monitor_iter;
  p | method_check_include_ghosts;

  p | method_gravity_type_super;

  p | method_inference_level_base;
  p | method_inference_level_array;
  p | method_inference_level_infer;
  p | method_inference_field_group;
  p | method_inference_overdensity_threshold;

  PUParray(p,initial_accretion_test_sink_position,3);
  PUParray(p,initial_accretion_test_sink_velocity,3);
  p | initial_accretion_test_sink_mass;
  p | initial_accretion_test_gas_density;
  p | initial_accretion_test_gas_pressure;
  p | initial_accretion_test_gas_radial_velocity;

  PUParray(p,initial_bb_test_center,3);
  PUParray(p,initial_bb_test_drift_velocity,3);
  p | initial_bb_test_mean_density;
  p | initial_bb_test_fluctuation_amplitude;
  p | initial_bb_test_truncation_radius;
  p | initial_bb_test_nominal_sound_speed;
  p | initial_bb_test_angular_rotation_velocity;
  p | initial_bb_test_external_density;

  p | method_turbulence_edot;

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
  read_initial_collapse_(p);
  read_initial_cosmology_(p);
  read_initial_hdf5_(p);
  read_initial_isolated_galaxy_(p);
  read_initial_merge_sinks_test_(p);
  read_initial_music_(p);
  read_initial_pm_(p);
  read_initial_sedov_(p);
  read_initial_sedov_random_(p);
  read_initial_turbulence_(p);

  // it's important for read_physics_
  read_physics_(p);

  // Method [sorted]

  read_method_check_(p);
  read_method_gravity_(p);
  read_method_inference_(p);
  read_method_turbulence_(p);

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

void EnzoConfig::read_initial_hdf5_(Parameters * p)
{
  const std::string name_initial = "Initial:hdf5:";

  initial_hdf5_max_level = p->value_integer (name_initial + "max_level", 0);
  initial_hdf5_format    = p->value_string  (name_initial + "format", "music");

  // Ensure hdf5 max level agrees with adapt max initial level.
  int adapt_max_level = p->value_integer("Adapt:max_level", 0);
  int adapt_max_initial_level = p->value_integer("Adapt:max_initial_level", adapt_max_level);
  if (initial_hdf5_max_level > 0 && adapt_max_initial_level != initial_hdf5_max_level) {
    ERROR2("Config::read_initial_hdf5_()",
    "The hdf5 max level (%d) should equal Adapt:max_initial_level (%d)",
    initial_hdf5_max_level,
    adapt_max_initial_level);
  }

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
    const int level           = p->value_integer (file_id + "level", 0);

    if (type == "particle") {

      const std::string attribute = p->value_string (file_id+"attribute","");

      initial_hdf5_particle_files.     push_back(file);
      initial_hdf5_particle_datasets.  push_back(dataset);
      initial_hdf5_particle_coords.    push_back(coords);
      initial_hdf5_particle_types.     push_back(name);
      initial_hdf5_particle_attributes.push_back(attribute);
      initial_hdf5_particle_levels.    push_back(level);

    } else if (type == "field") {

      initial_hdf5_field_files.        push_back(file);
      initial_hdf5_field_datasets.     push_back(dataset);
      initial_hdf5_field_names.        push_back(name);
      initial_hdf5_field_coords.       push_back(coords);
      initial_hdf5_field_levels.       push_back(level);

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

void EnzoConfig::read_initial_bcenter_(Parameters * p)
{
  // VL+CT b-field initialization
  initial_bcenter_update_etot = p->value_logical
    ("Initial:vlct_bfield:update_etot",false);
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
  method_check_monitor_iter   = p->value_integer("monitor_iter",0);
  method_check_include_ghosts = p->value_logical("include_ghosts",false);
}

//----------------------------------------------------------------------

void EnzoConfig::read_method_gravity_(Parameters * p)
{
  method_gravity_type_super = p->value_string
    ("Method:gravity:type_super","accelerations");
}

//-------------------------------------------------------------

void EnzoConfig::read_method_inference_(Parameters* p)
{
  p->group_set(0,"Method");
  p->group_push("inference");

  method_inference_level_base = p->value_integer ("level_base");
  method_inference_level_array = p->value_integer ("level_array");
  method_inference_level_infer = p->value_integer ("level_infer");

  const int rank = p->value_integer("Mesh:root_rank",0);

  method_inference_field_group = p->value_string  ("field_group");

  method_inference_overdensity_threshold = p->value_float
    ("Method:inference:overdensity_threshold",0.0);
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
               "there's no support for building an eos of type \"%s\".",
               type.c_str());
      }
    }


    if (legacy_gamma_specified) {
      double gamma = p->value_float("Field:gamma", -1.0);
      if (gamma <= 1.0) {
        std::string isothermal_name = EnzoEOSIsothermal::name();
        ERROR2("parse_eos_choice_",
               "\"Field:gamma\" is a legacy parameter that will be removed. "
               "It has an invalid value of 1 or smaller. If you want to "
               "initialize an \"%s\" EOS, you should delete this parameter & "
               "assign Physics:fluid_props:eos:type a value of \"%s\"",
               isothermal_name.c_str(), isothermal_name.c_str());
      }
      WARNING1("parse_eos_choice_",
               "\"Field:gamma\" is a legacy parameter that will be removed. "
               "It is being used to configure an \"%s\" EOS. Going forward, "
               "set parameters in the \"Physics:fluid_props:eos\" parameter "
               "group instead.",
               ideal_name.c_str());
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
