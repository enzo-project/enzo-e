// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Config.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-02
/// @brief    [\ref Parameters] Declaration of the Config class
///
/// Last review of parameters was 2015-09-10 with revision -r 3836 

#ifndef PARAMETERS_CONFIG_HPP
#define PARAMETERS_CONFIG_HPP

/* Maximum number of fields in any field list in the configuration file */
#define MAX_FIELDS      30

/* Maximum number of method groups in the configuration file */
#define MAX_METHOD_GROUPS  10

/* Maximum number of output file groups specified in the configuration file */
#define MAX_OUTPUT_GROUPS 20

/* Maximum number of schedules */
#define MAX_SCHEDULE 30

/* Maximum number of Boundary groups specified in the configuration file */

#define MAX_BOUNDARY 50

class Parameters;

class Config : public PUP::able {

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

  // Adapt

  int                        num_adapt;
  std::vector <std::string>  adapt_list;
  int                        adapt_interval;
  int                        adapt_min_face_rank;
  std::vector <std::string>  adapt_type;
  std::vector 
  < std::vector<std::string> > adapt_field_list;
  std::vector <double>       adapt_min_refine;
  std::vector <double>       adapt_max_coarsen;
  std::vector <double>       adapt_min_refine2;
  std::vector <double>       adapt_max_coarsen2;
  std::vector <int>          adapt_max_level;
  std::vector <double>       adapt_level_exponent;
  std::vector <char>         adapt_include_ghosts;
  std::vector <std::string>  adapt_output;
  
  // Balance (dynamic load balancing)

  int                        balance_schedule_index;

  // Boundary

  int                        num_boundary;
  std::string                boundary_list[MAX_BOUNDARY];
  std::string                boundary_type[MAX_BOUNDARY];
  int                        boundary_axis[MAX_BOUNDARY];
  int                        boundary_face[MAX_BOUNDARY];
  int                        boundary_mask[MAX_BOUNDARY];
  std::vector<std::string>   boundary_field_list[MAX_BOUNDARY];

  // Domain

  double                     domain_lower[3];
  double                     domain_upper[3];

  // Field

  int                        num_fields;
  std::vector<std::string>   field_list;
  std::map<std::string,int>  field_index;
  int                        field_alignment;
  std::vector<int>           field_centering [3];
  int                        field_ghost_depth[3];
  int                        field_padding;
  int                        field_precision;
  std::string                field_prolong;
  std::string                field_restrict;
  std::vector< std::vector<std::string> >  field_group_list;

  // Initial

  int                        num_initial;
  std::vector<std::string>   initial_list;
  int                        initial_cycle;
  double                     initial_time;

  // Memory

  bool                       memory_active;

  // Mesh

  int                        mesh_root_blocks[3];
  int                        mesh_root_rank;
  int                        mesh_root_size[3];
  int                        mesh_min_level;
  int                        mesh_max_level;

  // Method

  int                        num_method;
  std::vector<std::string>   method_list;
  int                        method_schedule_index [MAX_METHOD_GROUPS];
  std::vector<double>        method_courant;

  // Monitor

  bool                       monitor_debug;
  bool                       monitor_verbose;

  // Output

  int                        num_output;
  std::vector<std::string>   output_list;
  std::string                output_type                 [MAX_OUTPUT_GROUPS];
  std::string                output_axis           [MAX_OUTPUT_GROUPS];
  int                        output_image_block_size     [MAX_OUTPUT_GROUPS];
  std::vector<double>        output_colormap       [MAX_OUTPUT_GROUPS];
  std::string                output_image_type           [MAX_OUTPUT_GROUPS];
  bool                       output_image_log            [MAX_OUTPUT_GROUPS];
  std::string                output_image_mesh_color     [MAX_OUTPUT_GROUPS];
  std::vector<int>           output_image_size           [MAX_OUTPUT_GROUPS];
  std::string                output_image_reduce_type    [MAX_OUTPUT_GROUPS];
  bool                       output_image_ghost          [MAX_OUTPUT_GROUPS];
  int                        output_image_face_rank      [MAX_OUTPUT_GROUPS];
  bool                       output_image_specify_bounds [MAX_OUTPUT_GROUPS];
  double                     output_image_min            [MAX_OUTPUT_GROUPS];
  double                     output_image_max            [MAX_OUTPUT_GROUPS];
  int                        output_schedule_index [MAX_OUTPUT_GROUPS];
  std::vector<std::string>   output_dir            [MAX_OUTPUT_GROUPS];
  int                        output_stride         [MAX_OUTPUT_GROUPS];
  std::vector<std::string>   output_field_list     [MAX_OUTPUT_GROUPS];
  std::vector<std::string>   output_particle_list     [MAX_OUTPUT_GROUPS];
  std::vector<std::string>   output_name           [MAX_OUTPUT_GROUPS];

  int                        index_schedule_;
  std::vector<double>        schedule_list  [MAX_SCHEDULE];
  std::string                schedule_type  [MAX_SCHEDULE];
  std::string                schedule_var   [MAX_SCHEDULE];
  double                     schedule_start [MAX_SCHEDULE];
  double                     schedule_stop  [MAX_SCHEDULE];
  double                     schedule_step  [MAX_SCHEDULE];

  // Particles

  int                        num_particles;  // number of particle types
  std::vector<std::string>   particle_list;
  std::map<std::string,int>  particle_index;
  std::vector<char>          particle_interleaved;
  std::vector< std::vector <std::string> > particle_attribute_name;
  std::vector< std::vector <std::string> > particle_attribute_type;
  int                        particle_batch_size;
  std::vector< std::vector<std::string> >  particle_group_list;

  // Performance

  std::vector<std::string>   performance_papi_counters;
  std::string                performance_name;
  int                        performance_stride;
  bool                       performance_warnings;

  // Restart

  std::string                restart_file;

  // Stopping

  int                        stopping_cycle;
  double                     stopping_time;
  double                     stopping_seconds;
  int                        stopping_interval;

  // Testing

  int                        testing_cycle_final;
  double                     testing_time_final;
  double                     testing_time_tolerance;

protected: // functions

  /// Read boundary-related values from the Parameters object
  void read_adapt_       ( Parameters * ) throw();
  void read_balance_     ( Parameters * ) throw();
  void read_boundary_    ( Parameters * ) throw();
  void read_domain_      ( Parameters * ) throw();
  void read_field_       ( Parameters * ) throw();
  void read_initial_     ( Parameters * ) throw();
  void read_memory_      ( Parameters * ) throw();
  void read_mesh_        ( Parameters * ) throw();
  void read_method_      ( Parameters * ) throw();
  void read_monitor_     ( Parameters * ) throw();
  void read_output_      ( Parameters * ) throw();
  void read_particle_    ( Parameters * ) throw();
  void read_performance_ ( Parameters * ) throw();
  void read_restart_     ( Parameters * ) throw();
  void read_stopping_    ( Parameters * ) throw();
  void read_testing_     ( Parameters * ) throw();

  int read_schedule_( Parameters * ,
		      const std::string group   );

};

#endif /* PARAMETERS_CONFIG_HPP */

