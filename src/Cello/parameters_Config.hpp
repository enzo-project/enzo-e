// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Config.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-02
/// @brief    [\ref Parameters] Declaration of the Config class
///

#ifndef PARAMETERS_CONFIG_HPP
#define PARAMETERS_CONFIG_HPP

#define MAX_FIELDS 10
#define MAX_FILE_GROUPS 10

class Parameters;

class Config {

  /// @class    Config
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] 

public: // interface

  /// Constructor
  Config() throw();

  /// Destructor
  ~Config() throw();

  /// Copy constructor
  Config(const Config & config) throw();

  /// Assignment operator
  Config & operator= (const Config & config) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
#endif

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();

public: // attributes

  // NOTE: change pup() function whenever attributes change

  std::string                boundary_type;

  double                     domain_lower[3];
  double                     domain_upper[3];

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

  int                        field_alignment;
  std::vector<bool>          field_centering [MAX_FIELDS];
  double                     field_courant;
  std::vector<std::string>   field_fields;
  int                        field_ghosts[3];;
  int                        field_padding;
  std::string                field_precision;
  bool                       field_refresh_corners;
  bool                       field_refresh_edges;

  int                        initial_cycle;
  std::string                initial_type;
  double                     initial_time;
  std::vector<std::string>   initial_name;
  std::vector<std::string>   initial_value [MAX_FIELDS];

  int                        mesh_root_blocks[3];
  int                        mesh_root_rank;
  int                        mesh_root_size[3];

  std::vector<std::string>   method_sequence;

  bool                       monitor_debug;

  std::vector<std::string>   output_file_groups;
  std::string                output_axis           [MAX_FILE_GROUPS];
  std::vector<double>        output_colormap       [MAX_FILE_GROUPS];
  std::vector<double>        output_colormap_alpha [MAX_FILE_GROUPS];
  std::vector<std::string>   output_field_list     [MAX_FILE_GROUPS];
  std::vector<std::string>   output_name           [MAX_FILE_GROUPS];
  std::vector<std::string>   output_dir            [MAX_FILE_GROUPS];
  std::vector<std::string>   output_schedule       [MAX_FILE_GROUPS];
  std::string                output_type           [MAX_FILE_GROUPS];

  bool                       physics_cosmology;
  double                     physics_cosmology_comoving_box_size;
  double                     physics_cosmology_hubble_constant_now;
  double                     physics_cosmology_initial_redshift;
  double                     physics_cosmology_max_expansion_rate;
  double                     physics_cosmology_omega_lamda_now;
  double                     physics_cosmology_omega_matter_now;
  double                     physics_gamma;

  int                        stopping_cycle;
  double                     stopping_time;

  int                        testing_cycle_final;
  double                     testing_time_final;

  std::string                timestep_type;

};

#endif /* PARAMETERS_CONFIG_HPP */

