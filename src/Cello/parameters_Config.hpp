// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Config.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-02
/// @brief    [\ref Parameters] Declaration of the Config class
///

#ifndef PARAMETERS_CONFIG_HPP
#define PARAMETERS_CONFIG_HPP

/* Maximum number of fields in any field list in the configuration file */
#define MAX_FIELDS      30

/* Maximum number of output file groups specified in the configuration file */
#define MAX_FILE_GROUPS 10

class Parameters;

class Config : public PUP::able 
{

  /// @class    Config
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] 

public: // interface

  /// empty constructor for charm++ pup()
  Config() throw() {}

  /// CHARM++ PUP::able declaration
  PUPable_decl(Config);

  /// CHARM++ migration constructor for PUP::able

  Config (CkMigrateMessage *m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read values from the Parameters object
  void read (Parameters * parameters) throw();

public: // attributes

  // NOTE: change pup() function whenever attributes change

  std::string                boundary_type;

  double                     domain_lower[3];
  double                     domain_upper[3];

  int                        num_fields;
  int                        field_alignment;
  std::vector<int>           field_centering [3];
  double                     field_courant;
  std::vector<std::string>   field_fields;
  int                        field_ghosts[3];
  int                        field_padding;
  int                        field_precision;
  int                        field_refresh_rank;

  int                        initial_cycle;
  std::string                initial_type;
  double                     initial_time;
  std::vector<std::string>   initial_name;
  // std::vector<std::string>   initial_value [MAX_FIELDS];
  int                        initial_max_level;

  bool                       memory_active;

  int                        mesh_root_blocks[3];
  int                        mesh_root_rank;
  int                        mesh_root_size[3];
  int                        mesh_max_level;
  bool                       mesh_balance;

  std::vector<std::string>   mesh_adapt_type ;
  std::vector<std::string>   mesh_adapt_fields ;
  double                     mesh_adapt_slope_min_refine;
  double                     mesh_adapt_slope_max_coarsen;
  double                     mesh_adapt_mass_min;
  double                     mesh_adapt_mass_min_overdensity;
  double                     mesh_adapt_mass_level_exponent;
  bool                       mesh_adapt_balance;
  int                        mesh_adapt_interval;
  std::string                mesh_sync_adapt_enter;
  std::string                mesh_sync_adapt_called;
  std::string                mesh_sync_adapt_next;
  std::string                mesh_sync_adapt_exit;
  std::string                mesh_sync_refresh_enter;
  std::string                mesh_sync_refresh_exit;

  std::vector<std::string>   method_sequence;

  bool                       monitor_debug;

  int                        num_file_groups;
  std::vector<std::string>   output_file_groups;
  std::string                output_type           [MAX_FILE_GROUPS];

  std::string                output_image_axis           [MAX_FILE_GROUPS];
  int                        output_image_block_size     [MAX_FILE_GROUPS];
  std::vector<double>        output_image_colormap_alpha [MAX_FILE_GROUPS];
  std::vector<double>        output_image_colormap       [MAX_FILE_GROUPS];
  std::string                output_image_type           [MAX_FILE_GROUPS];
  bool                       output_image_log            [MAX_FILE_GROUPS];
  std::string                output_image_mesh_color     [MAX_FILE_GROUPS];
  std::vector<int>           output_image_size           [MAX_FILE_GROUPS];
  std::string                output_image_reduce_type    [MAX_FILE_GROUPS];
  bool                       output_image_ghost          [MAX_FILE_GROUPS];
  int                        output_image_face_rank      [MAX_FILE_GROUPS];
  bool                       output_image_specify_bounds [MAX_FILE_GROUPS];
  double                     output_image_min            [MAX_FILE_GROUPS];
  double                     output_image_max            [MAX_FILE_GROUPS];
  std::vector<std::string>   output_dir            [MAX_FILE_GROUPS];
  int                        output_stride         [MAX_FILE_GROUPS];
  std::vector<std::string>   output_field_list     [MAX_FILE_GROUPS];
  std::vector<std::string>   output_name           [MAX_FILE_GROUPS];
  std::string                output_schedule_type  [MAX_FILE_GROUPS];
  std::string                output_schedule_var   [MAX_FILE_GROUPS];
  double                     output_schedule_start [MAX_FILE_GROUPS];
  double                     output_schedule_stop  [MAX_FILE_GROUPS];
  double                     output_schedule_step  [MAX_FILE_GROUPS];
  std::vector<double>        output_schedule_list  [MAX_FILE_GROUPS];

  std::vector<std::string>   performance_papi_counters;
  std::string                performance_name;
  int                        performance_stride;
  bool                       performance_warnings;

  std::string                prolong_type;
  std::string                restrict_type;

  int                        stopping_cycle;
  double                     stopping_time;
  int                        stopping_interval;

  int                        testing_cycle_final;
  double                     testing_time_final;

  std::string                timestep_type;

};

#endif /* PARAMETERS_CONFIG_HPP */

