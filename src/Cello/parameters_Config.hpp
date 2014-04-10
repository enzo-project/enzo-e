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

/* Maximum number of schedules */
#define MAX_SCHEDULE 10

/* Maximum number of adapt groups specified in the configuration file */
#define MAX_ADAPT 10

/* Maximum number of Boundary groups specified in the configuration file */

#define MAX_BOUNDARY 50

class Parameters;

class Config : public PUP::able 
{

  /// @class    Config
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] 

public: // interface

  /// empty constructor for charm++ pup()
  Config() throw() 
  : index_schedule_(0) {}

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

  int                        num_boundary;
  std::string                boundary_list[MAX_BOUNDARY];
  std::string                boundary_type[MAX_BOUNDARY];
  int                        boundary_axis[MAX_BOUNDARY];
  int                        boundary_face[MAX_BOUNDARY];
  int                        boundary_mask[MAX_BOUNDARY];
  std::vector<std::string>   boundary_field_list[MAX_BOUNDARY];

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
  int                        initial_max_level;

  bool                       memory_active;

  int                        mesh_root_blocks[3];
  int                        mesh_root_rank;
  int                        mesh_root_size[3];
  int                        mesh_max_level;

  int                        adapt_interval;
  bool                       adapt_balance;
  int                        num_adapt;
  std::string                adapt_list[MAX_ADAPT];
  std::string                adapt_type[MAX_ADAPT];
  std::vector<std::string>   adapt_field_list[MAX_ADAPT];
  double                     adapt_min_refine[MAX_ADAPT];
  double                     adapt_max_coarsen[MAX_ADAPT];
  double                     adapt_level_exponent[MAX_ADAPT];

  std::string                control_sync_adapt_enter;
  std::string                control_sync_adapt_called;
  std::string                control_sync_adapt_end;
  std::string                control_sync_adapt_next;
  std::string                control_sync_adapt_exit;
  std::string                control_sync_compute_enter;
  std::string                control_sync_compute_exit;
  std::string                control_sync_exit;
  std::string                control_sync_output_enter;
  std::string                control_sync_output_exit;
  std::string                control_sync_refresh_enter;
  std::string                control_sync_refresh_exit;
  std::string                control_sync_stopping_enter;
  std::string                control_sync_stopping_exit;

  std::vector<std::string>   method_sequence;

  bool                       monitor_debug;

  int                        num_file_groups;
  std::vector<std::string>   output_file_groups;
  std::string                output_type           [MAX_FILE_GROUPS];

  std::string                output_image_axis           [MAX_FILE_GROUPS];
  int                        output_image_block_size     [MAX_FILE_GROUPS];
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
  int                        output_schedule_index [MAX_FILE_GROUPS];
  std::vector<std::string>   output_dir            [MAX_FILE_GROUPS];
  int                        output_stride         [MAX_FILE_GROUPS];
  std::vector<std::string>   output_field_list     [MAX_FILE_GROUPS];
  std::vector<std::string>   output_name           [MAX_FILE_GROUPS];
  int                        index_schedule_;
  std::string                output_schedule_type  [MAX_SCHEDULE];
  std::string                output_schedule_var   [MAX_SCHEDULE];
  double                     output_schedule_start [MAX_SCHEDULE];
  double                     output_schedule_stop  [MAX_SCHEDULE];
  double                     output_schedule_step  [MAX_SCHEDULE];
  std::vector<double>        output_schedule_list  [MAX_SCHEDULE];

  std::vector<std::string>   performance_papi_counters;
  std::string                performance_name;
  int                        performance_stride;
  bool                       performance_warnings;

  // schedule for turning on / off Projections monitoring

  std::string                projections_schedule_on_type;
  std::string                projections_schedule_on_var;
  double                     projections_schedule_on_start;
  double                     projections_schedule_on_stop;
  double                     projections_schedule_on_step;
  std::vector<double>        projections_schedule_on_list;
  std::string                projections_schedule_off_type;
  std::string                projections_schedule_off_var;
  double                     projections_schedule_off_start;
  double                     projections_schedule_off_stop;
  double                     projections_schedule_off_step;
  std::vector<double>        projections_schedule_off_list;

  std::string                prolong_type;
  std::string                restrict_type;

  int                        stopping_cycle;
  double                     stopping_time;
  int                        stopping_interval;

  int                        testing_cycle_final;
  double                     testing_time_final;

  int                        num_timestep;
  std::vector<std::string>   timestep_type;

protected: // functions

  /// Read boundary-related values from the Parameters object
  void read_boundary_    (Parameters * parameters) throw();
  void read_domain_      (Parameters * parameters) throw();
  void read_field_       (Parameters * parameters) throw();
  void read_initial_     (Parameters * parameters) throw();
  void read_memory_      (Parameters * parameters) throw();
  void read_mesh_        (Parameters * parameters) throw();
  void read_control_     (Parameters * parameters) throw();
  void read_adapt_       (Parameters * parameters) throw();
  void read_method_      (Parameters * parameters) throw();
  void read_monitor_     (Parameters * parameters) throw();
  void read_output_      (Parameters * parameters) throw();
  void read_performance_ (Parameters * parameters) throw();
  void read_stopping_    (Parameters * parameters) throw();
  void read_testing_     (Parameters * parameters) throw();
  void read_timestep_    (Parameters * parameters) throw();

  void read_schedule_(Parameters * parameters,
		      const std::string group,
		      std::string * type,
		      std::string * var,
		      double * start,
		      double * stop,
		      double * step,
		      std::vector<double> & list);

  int read_schedule_new_(Parameters * parameters,
		      const std::string group   );

};

#endif /* PARAMETERS_CONFIG_HPP */

