// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Config.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-02
/// @brief    [\ref Parameters] Declaration of the Config class
///

#ifndef PARAMETERS_CONFIG_HPP
#define PARAMETERS_CONFIG_HPP

#define MAX_FIELDS      30
#define MAX_FILE_GROUPS 10

class Parameters;

class Config {

  /// @class    Config
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] 

public: // interface

  // /// Constructor
  // Config() throw();

  // /// Destructor
  // ~Config() throw();

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

  int                        num_fields;
  int                        field_alignment;
  std::vector<int>           field_centering [3];
  double                     field_courant;
  std::vector<std::string>   field_fields;
  int                        field_ghosts[3];;
  int                        field_padding;
  int                        field_precision;
  bool                       field_refresh_corners;
  bool                       field_refresh_edges;
  bool                       field_refresh_faces;

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

  int                        num_file_groups;
  std::vector<std::string>   output_file_groups;
  std::string                output_type           [MAX_FILE_GROUPS];

  std::string                output_image_axis           [MAX_FILE_GROUPS];
  std::vector<double>        output_image_colormap_alpha [MAX_FILE_GROUPS];
  std::vector<double>        output_image_colormap       [MAX_FILE_GROUPS];
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

  int                        stopping_cycle;
  double                     stopping_time;

  int                        testing_cycle_final;
  double                     testing_time_final;

  std::string                timestep_type;

};

#endif /* PARAMETERS_CONFIG_HPP */

