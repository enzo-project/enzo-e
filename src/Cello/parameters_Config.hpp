// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Config.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-02
/// @brief    [\ref Parameters] Declaration of the Config class
///
/// Last review of parameters was 2015-09-10 with revision -r 3836 

#ifndef PARAMETERS_CONFIG_HPP
#define PARAMETERS_CONFIG_HPP

class Parameters;

class Config : public PUP::able {

  /// @class    Config
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] 

public: // interface

  /// empty constructor for charm++ pup()
  Config() throw() 
    : PUP::able(),
      num_adapt(0),
      adapt_list(),
      adapt_interval(0),
      adapt_min_face_rank(0),
      adapt_type(),
      adapt_field_list(),
      adapt_min_refine(),
      adapt_max_coarsen(),
      adapt_min_refine2(),
      adapt_max_coarsen2(),
      adapt_max_level(),
      adapt_level_exponent(),
      adapt_include_ghosts(),
      adapt_output(),
      balance_schedule_index(0),
      num_boundary(0),
      boundary_list(),
      boundary_type(),
      boundary_axis(),
      boundary_face(),
      boundary_mask(),
      boundary_field_list(),
      num_fields(0),
      field_list(),
      field_index(),
      field_alignment(0),
      field_padding(0),
      field_precision(0),
      field_prolong(""),
      field_restrict(""),
      field_group_list(),
      num_initial(0),
      initial_list(),
      initial_cycle(0),
      initial_time(0.0),
      initial_trace_field(""),
      initial_trace_mpp(0.0),
      initial_trace_dx(0),
      initial_trace_dy(0),
      initial_trace_dz(0),
      memory_active(false),
      memory_warning_mb(0.0),
      memory_limit_gb(0.0),
      mesh_root_rank(0),
      mesh_min_level(0),
      mesh_max_level(0),
      num_method(0),
      method_courant_global(1.0),
      method_list(),
      method_schedule_index(),
      method_courant(),
      method_timestep(),
      monitor_debug(false),
      monitor_verbose(false),
      num_output(0),
      output_list(),
      output_type(),
      output_axis(),
      output_image_block_size(),
      output_image_lower(),
      output_image_upper(),
      output_colormap(),
      output_image_type(),
      output_image_log(),
      output_image_abs(),
      output_image_mesh_color(),
      output_image_color_particle_attribute(),
      output_image_size(),
      output_image_reduce_type(),
      output_image_ghost(),
      output_image_face_rank(),
      output_image_min(),
      output_image_max(),
      output_schedule_index(),
      output_max_level(),
      output_min_level(),
      output_leaf_only(),
      output_dir(),
      output_stride(),
      output_field_list(),
      output_particle_list(),
      output_name(),
      index_schedule_(0),
      schedule_list(),
      schedule_type(),
      schedule_var(),
      schedule_start(),
      schedule_stop(),
      schedule_step(),
      num_particles(0),
      particle_list(),
      particle_index(),
      particle_interleaved(),
      particle_constant_name(),
      particle_constant_type(),
      particle_constant_value(),
      particle_attribute_name(),
      particle_attribute_type(),
      particle_batch_size(0),
      particle_group_list(),
      performance_papi_counters(),
      performance_warnings(false),
      restart_file(""),
      num_solvers(),
      solver_list(),
      solver_index(),
      solver_type(),
      solver_iter_max(),
      solver_res_tol(),
      solver_diag_precon(),
      solver_monitor_iter(),
      solver_restrict(),
      solver_prolong(),
      solver_min_level(),
      solver_max_level(),
      stopping_cycle(0),
      stopping_time(0.0),
      stopping_seconds(0.0),
      stopping_interval(0),
      testing_cycle_final(0),
      testing_time_final(0.0),
      testing_time_tolerance(0.0)
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(Config);

  /// CHARM++ migration constructor for PUP::able

  Config (CkMigrateMessage *m)
    : PUP::able(m),
      num_adapt(0),
      adapt_list(),
      adapt_interval(0),
      adapt_min_face_rank(0),
      adapt_type(),
      adapt_field_list(),
      adapt_min_refine(),
      adapt_max_coarsen(),
      adapt_min_refine2(),
      adapt_max_coarsen2(),
      adapt_max_level(),
      adapt_level_exponent(),
      adapt_include_ghosts(),
      adapt_output(),
      balance_schedule_index(0),
      num_boundary(0),
      boundary_list(),
      boundary_type(),
      boundary_axis(),
      boundary_face(),
      boundary_mask(),
      boundary_field_list(),
      num_fields(0),
      field_list(),
      field_index(),
      field_alignment(0),
      field_padding(0),
      field_precision(0),
      field_prolong(""),
      field_restrict(""),
      field_group_list(),
      num_initial(0),
      initial_list(),
      initial_cycle(0),
      initial_time(0.0),
      initial_trace_field(""),
      initial_trace_mpp(0.0),
      initial_trace_dx(0),
      initial_trace_dy(0),
      initial_trace_dz(0),
      memory_active(false),
      memory_warning_mb(0.0),
      memory_limit_gb(0.0),
      mesh_root_rank(0),
      mesh_min_level(0),
      mesh_max_level(0),
      num_method(0),
      method_courant_global(1.0),
      method_list(),
      method_schedule_index(),
      method_courant(),
      method_timestep(),
      monitor_debug(false),
      monitor_verbose(false),
      num_output(0),
      output_list(),
      output_type(),
      output_axis(),
      output_image_block_size(),
      output_image_lower(),
      output_image_upper(),
      output_colormap(),
      output_image_type(),
      output_image_log(),
      output_image_abs(),
      output_image_mesh_color(),
      output_image_color_particle_attribute(),
      output_image_size(),
      output_image_reduce_type(),
      output_image_ghost(),
      output_image_face_rank(),
      output_image_min(),
      output_image_max(),
      output_schedule_index(),
      output_max_level(),
      output_min_level(),
      output_leaf_only(),
      output_dir(),
      output_stride(),
      output_field_list(),
      output_particle_list(),
      output_name(),
      index_schedule_(0),
      schedule_list(),
      schedule_type(),
      schedule_var(),
      schedule_start(),
      schedule_stop(),
      schedule_step(),
      num_particles(0),
      particle_list(),
      particle_index(),
      particle_interleaved(),
      particle_constant_name(),
      particle_constant_type(),
      particle_constant_value(),
      particle_attribute_name(),
      particle_attribute_type(),
      particle_batch_size(0),
      particle_group_list(),
      performance_papi_counters(),
      performance_warnings(false),
      restart_file(""),
      num_solvers(),
      solver_list(),
      solver_index(),
      solver_type(),
      solver_iter_max(),
      solver_res_tol(),
      solver_diag_precon(),
      solver_monitor_iter(),
      solver_restrict(),
      solver_prolong(),
      solver_min_level(),
      solver_max_level(),
      stopping_cycle(0),
      stopping_time(0.0),
      stopping_seconds(0.0),
      stopping_interval(0),
      testing_cycle_final(0),
      testing_time_final(0.0),
      testing_time_tolerance(0.0)
  {
    for (int axis=0; axis<3; axis++) {
      domain_lower[axis] = 0.0;
      domain_upper[axis] = 0.0;
      field_ghost_depth[axis] = 0;
      mesh_root_blocks[axis] = 0;
      mesh_root_size[axis] = 0;
    }
  }

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
  std::vector<std::string>   boundary_list;
  std::vector<std::string>   boundary_type;
  std::vector<int>           boundary_axis;
  std::vector<int>           boundary_face;
  std::vector<int>           boundary_mask;
  std::vector<std::vector<std::string> >  
                             boundary_field_list;

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

  std::string                initial_trace_field;
  double                     initial_trace_mpp;
  int                        initial_trace_dx;
  int                        initial_trace_dy;
  int                        initial_trace_dz;

  // Memory

  bool                       memory_active;
  double                     memory_warning_mb;
  double                     memory_limit_gb;

  // Mesh

  int                        mesh_root_blocks[3];
  int                        mesh_root_rank;
  int                        mesh_root_size[3];
  int                        mesh_min_level;
  int                        mesh_max_level;

  // Method

  int                        num_method;
  double                     method_courant_global;
  std::vector<std::string>   method_list;
  std::vector<int>           method_schedule_index;
  std::vector<double>        method_courant;
  std::vector<double>        method_timestep;

  // Monitor

  bool                       monitor_debug;
  bool                       monitor_verbose;

  // Output

  int                         num_output;
  std::vector <std::string>   output_list;
  std::vector < std::string > output_type;
  std::vector < std::string > output_axis;
  std::vector < int >         output_image_block_size;
  std::vector < std::vector <double> > output_image_lower;
  std::vector < std::vector <double> > output_image_upper;
  std::vector < std::vector <double> > output_colormap;
  std::vector < std::string > output_image_type;
  std::vector < char >        output_image_log;
  std::vector < char >        output_image_abs;
  std::vector < std::string > output_image_mesh_color;
  std::vector < std::string > output_image_color_particle_attribute;
  std::vector < std::vector <int> > output_image_size;
  std::vector < std::string>  output_image_reduce_type;
  std::vector < char>         output_image_ghost;
  std::vector < int >         output_image_face_rank;
  std::vector < double>       output_image_min;
  std::vector < double>       output_image_max;
  std::vector < int >         output_schedule_index;
  std::vector < int >         output_max_level;
  std::vector < int >         output_min_level;
  std::vector < char >        output_leaf_only;
  std::vector < std::vector <std::string> >  output_dir;
  std::vector < int >         output_stride;
  std::vector < std::vector <std::string> >  output_field_list;
  std::vector < std::vector <std::string> > output_particle_list;
  std::vector < std::vector <std::string> >  output_name;
  int                        index_schedule_;
  std::vector< std::vector<double> > schedule_list;
  std::vector< std::string > schedule_type;
  std::vector< std::string > schedule_var;
  std::vector< double >      schedule_start;
  std::vector< double >      schedule_stop;
  std::vector< double >      schedule_step;

  // Particles

  int                        num_particles;  // number of particle types
  std::vector<std::string>   particle_list;
  std::map<std::string,int>  particle_index;
  std::vector<char>          particle_interleaved;

  std::vector< std::vector <std::string> > particle_constant_name;
  std::vector< std::vector <std::string> > particle_constant_type;
  std::vector< std::vector <double> >      particle_constant_value;

  std::vector< std::vector <std::string> > particle_attribute_name;
  std::vector< std::vector <std::string> > particle_attribute_type;
  std::vector <int>          particle_attribute_position[3];
  std::vector <int>          particle_attribute_velocity[3];

  int                        particle_batch_size;
  std::vector< std::vector<std::string> >  particle_group_list;

  // Performance

  std::vector<std::string>   performance_papi_counters;
  bool                       performance_warnings;

  // Restart

  std::string                restart_file;

  // Solvers

  int                        num_solvers;
  std::vector<std::string>   solver_list;
  std::map<std::string,int>  solver_index;
  std::vector<std::string>   solver_type;
  std::vector<int>           solver_iter_max;
  std::vector<double>        solver_res_tol;
  std::vector<char>          solver_diag_precon;
  std::vector<int>           solver_monitor_iter;
  std::vector<std::string>   solver_restrict;
  std::vector<std::string>   solver_prolong;
  std::vector<int>           solver_min_level;
  std::vector<int>           solver_max_level;

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
  void read_solver_      ( Parameters * ) throw();
  void read_stopping_    ( Parameters * ) throw();
  void read_testing_     ( Parameters * ) throw();

  int read_schedule_( Parameters * ,
		      const std::string group   );

};

extern Config g_config;

#endif /* PARAMETERS_CONFIG_HPP */

