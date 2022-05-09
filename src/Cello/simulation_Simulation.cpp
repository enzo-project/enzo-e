// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "main.hpp"

#include "simulation.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_SIMULATION
// #define DEBUG_MSG_REFINE

Simulation::Simulation
(
 const char *   parameter_file,
 int            n
 )
/// Initialize the Simulation object
:
#if defined(CELLO_DEBUG) || defined(CELLO_VERBOSE)
  fp_debug_(NULL),
#endif
  factory_(NULL),
  parameters_(&g_parameters),
  parameter_file_(parameter_file),
  rank_(0),
  cycle_(0),
  cycle_watch_(-1),
  time_(0.0),
  dt_(0),
  stop_(false),
  phase_(phase_unknown),
  config_(&g_config),
  problem_(NULL),
  timer_(),
  performance_(NULL),
#ifdef CONFIG_USE_PROJECTIONS
  projections_tracing_(true),
  projections_schedule_on_(NULL),
  projections_schedule_off_(NULL),
#endif
  schedule_balance_(NULL),
  monitor_(NULL),
  hierarchy_(NULL),
  scalar_descr_long_double_(NULL),
  scalar_descr_double_(NULL),
  scalar_descr_int_(NULL),
  scalar_descr_long_long_(NULL),
  scalar_descr_sync_(NULL),
  scalar_descr_void_(NULL),
  scalar_descr_index_(NULL),
  field_descr_(NULL),
  particle_descr_(NULL),
  sync_output_begin_(),
  sync_output_write_(),
  sync_restart_created_(),
  sync_restart_next_(),
  refresh_list_(),
  index_output_(-1),
  num_solver_iter_(),
  max_solver_iter_(),
  restart_directory_(),
  restart_num_files_(),
  restart_stream_file_list_()
  
{
  for (int i=0; i<256; i++) dir_checkpoint_[i] = '\0';
#ifdef DEBUG_SIMULATION
  CkPrintf ("%d DEBUG_SIMULATION Simulation(parameter_file,n)\n",CkMyPe());
  fflush(stdout);
  char name[40];
  sprintf (name,"parameters-%02d.text",CkMyPe());
  parameters_->write(name);
#endif
  
  debug_open();

  monitor_ = Monitor::instance();
#ifdef CELLO_DEBUG
  monitor_->set_mode(monitor_mode_all);
#else
  monitor_->set_mode(monitor_mode_root);
#endif

}

//----------------------------------------------------------------------

Simulation::Simulation()
  :
#if defined(CELLO_DEBUG) || defined(CELLO_VERBOSE)
  fp_debug_(NULL),
#endif
  factory_(NULL),
  parameters_(&g_parameters),
  parameter_file_(""),
  rank_(0),
  cycle_(0),
  cycle_watch_(-1),
  time_(0.0),
  dt_(0),
  stop_(false),
  phase_(phase_unknown),
  config_(&g_config),
  problem_(NULL),
  timer_(),
  performance_(NULL),
#ifdef CONFIG_USE_PROJECTIONS
  projections_tracing_(true),
  projections_schedule_on_(NULL),
  projections_schedule_off_(NULL),
#endif
  schedule_balance_(NULL),
  monitor_(NULL),
  hierarchy_(NULL),
  scalar_descr_long_double_(NULL),
  scalar_descr_double_(NULL),
  scalar_descr_int_(NULL),
  scalar_descr_long_long_(NULL),
  scalar_descr_sync_(NULL),
  scalar_descr_void_(NULL),
  scalar_descr_index_(NULL),
  field_descr_(NULL),
  particle_descr_(NULL),
  sync_output_begin_(),
  sync_output_write_(),
  sync_restart_created_(),
  sync_restart_next_(),
  refresh_list_(),
  index_output_(-1),
  num_solver_iter_(),
  max_solver_iter_(),
  restart_directory_(),
  restart_num_files_(),
  restart_stream_file_list_()
{
  for (int i=0; i<256; i++) dir_checkpoint_[i] = '\0';
#ifdef DEBUG_SIMULATION
  CkPrintf ("%d DEBUG_SIMULATION Simulation()\n",CkMyPe());
  fflush(stdout);
#endif  
  TRACE("Simulation()");
}

//----------------------------------------------------------------------

Simulation::Simulation (CkMigrateMessage *m)
  : CBase_Simulation(m),
#if defined(CELLO_DEBUG) || defined(CELLO_VERBOSE)
    fp_debug_(NULL),
#endif
    factory_(NULL),
    parameters_(&g_parameters),
    parameter_file_(""),
    rank_(0),
    cycle_(0),
    cycle_watch_(-1),
    time_(0.0),
    dt_(0),
    stop_(false),
    phase_(phase_unknown),
    config_(&g_config),
    problem_(NULL),
    timer_(),
    performance_(NULL),
#ifdef CONFIG_USE_PROJECTIONS
    projections_tracing_(true),
    projections_schedule_on_(NULL),
    projections_schedule_off_(NULL),
#endif
    schedule_balance_(NULL),
    monitor_(NULL),
    hierarchy_(NULL),
    scalar_descr_long_double_(NULL),
    scalar_descr_double_(NULL),
    scalar_descr_int_(NULL),
    scalar_descr_long_long_(NULL),
    scalar_descr_sync_(NULL),
    scalar_descr_void_(NULL),
    scalar_descr_index_(NULL),
    field_descr_(NULL),
    particle_descr_(NULL),
    sync_output_begin_(),
    sync_output_write_(),
    sync_restart_created_(),
    sync_restart_next_(),
    refresh_list_(),
    index_output_(-1),
    num_solver_iter_(),
    max_solver_iter_(),
    restart_directory_(),
    restart_num_files_(),
    restart_stream_file_list_()

{
  for (int i=0; i<256; i++) dir_checkpoint_[i] = '\0';
#ifdef DEBUG_SIMULATION
  CkPrintf ("%d DEBUG_SIMULATION Simulation(msg)\n",CkMyPe());
  fflush(stdout);
#endif  
  TRACE("Simulation(CkMigrateMessage)");
}

//----------------------------------------------------------------------

Simulation::~Simulation()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Simulation::pup (PUP::er &p)
{
#ifdef DEBUG_SIMULATION
  CkPrintf ("%d DEBUG_SIMULATION Simulation::pup()\n",CkMyPe());
  fflush(stdout);
#endif  
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  CBase_Simulation::pup(p);

  bool up = p.isUnpacking();

  if (up) debug_open();

  p | factory_; // PUP::able

  p | config_;

  p | parameter_file_;

  p | rank_; 
  p | cycle_;
  p | cycle_watch_;
  p | time_;
  p | dt_;
  p | stop_;
  p | phase_;

  p | problem_; // PUPable

  if (up) performance_ = new Performance;
  p | *performance_;

  if (up) monitor_ = Monitor::instance();
  p | *monitor_;

  if (up) hierarchy_ = new Hierarchy;
  p | *hierarchy_;

  if (up) scalar_descr_long_double_ = new ScalarDescr;
  p | *scalar_descr_long_double_;
  if (up) scalar_descr_double_ = new ScalarDescr;
  p | *scalar_descr_double_;
  if (up) scalar_descr_int_ = new ScalarDescr;
  p | *scalar_descr_int_;
  if (up) scalar_descr_long_long_ = new ScalarDescr;
  p | *scalar_descr_long_long_;
  if (up) scalar_descr_sync_ = new ScalarDescr;
  p | *scalar_descr_sync_;
  if (up) scalar_descr_void_ = new ScalarDescr;
  p | *scalar_descr_void_;
  if (up) scalar_descr_index_ = new ScalarDescr;
  p | *scalar_descr_index_;

  if (up) field_descr_ = new FieldDescr;
  p | *field_descr_;

  if (up) particle_descr_ = new ParticleDescr;
  p | *particle_descr_;

  if (up && (phase_ == phase_restart)) {
    monitor_->header();
    monitor_->print ("Simulation","restarting");
  }

  p | sync_output_begin_;
  p | sync_output_write_;
  p | sync_restart_created_;
  p | sync_restart_next_;

  if (up) sync_output_begin_.set_stop(0);
  if (up) sync_output_write_.set_stop(0);

#ifdef CONFIG_USE_PROJECTIONS
  p | projections_tracing_;
  p | projections_schedule_on_;
  p | projections_schedule_off_;
#endif

  p | schedule_balance_;

  p | refresh_list_;
  p | refresh_name_;

  PUParray(p,dir_checkpoint_,256);

#ifdef BYPASS_CHARM_MEM_LEAK
  ASSERT1("Simulation::pup()",
	  "msg_refine_map_ is assumed to be empty but has size %lu",
	  msg_refine_map_.size(),
	  (msg_refine_map_.size() == 0));

  //  p | msg_refine_map_;
#endif
  p | index_output_;
  p | num_solver_iter_;
  p | max_solver_iter_;
  p | restart_directory_;
  p | restart_num_files_;

}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  TRACE0;

  performance_->stop_region(perf_simulation);

  performance_->end();

}

//----------------------------------------------------------------------

#ifdef BYPASS_CHARM_MEM_LEAK

void Simulation::p_get_msg_refine(Index index)
{
  MsgRefine * msg = get_msg_refine(index);
#ifdef DEBUG_MSG_REFINE  
  CkPrintf ("%d DEBUG_MSG_REFINE sending %p\n",CkMyPe(),msg);
#endif
  hierarchy_->block_array()[index].p_set_msg_refine(msg);
}

//----------------------------------------------------------------------

void Simulation::set_msg_refine(Index index, MsgRefine * msg)
{
  if (msg_refine_map_[index] != NULL) {
   
    int v3[3];
    index.values(v3);
    ASSERT3 ("Simulation::p_set_msg_refine",
	    "index %08x %08x %08x is already in the msg_refine mapping",
	    v3[0],v3[1],v3[2],
	    (msg == NULL));
  }
  msg_refine_map_[index] = msg;
}

//----------------------------------------------------------------------

MsgRefine * Simulation::get_msg_refine(Index index)
{
  int v3[3];
  index.values(v3);
  MsgRefine * msg = msg_refine_map_[index];
  if (msg == NULL) {
    int v3[3];
    index.values(v3);
    
    ASSERT3 ("Simulation::get_msg_refine",
	    "index %08x %08x %08x is not in the msg_refine mapping",
	    v3[0],v3[1],v3[2],
	    (msg != NULL));
  }
  msg_refine_map_.erase(index);
  return msg;
}

#endif
//======================================================================

void Simulation::initialize_simulation_() throw()
{

#ifdef DEBUG_SIMULATION
  CkPrintf ("%d DEBUG_SIMULATION Simulation::initialize_simulation_()\n",CkMyPe());
  fflush(stdout);
#endif  

  rank_ = config_->mesh_root_rank;

  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be specified",
	  rank_ != 0);
  
  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be 1, 2, or 3",
	  (1 <= rank_) && (rank_ <= 3));

  cycle_ = config_->initial_cycle;
  cycle_watch_ = cycle_ - 1;
  time_  = config_->initial_time;
  dt_ = 0;
}

//----------------------------------------------------------------------

void Simulation::initialize_memory_() throw()
{
  Memory * memory = Memory::instance();
  if (memory) {
    memory->set_active(config_->memory_active);
    memory->set_warning_mb (config_->memory_warning_mb);
    memory->set_limit_gb (config_->memory_limit_gb);
  }
  
}
//----------------------------------------------------------------------

void Simulation::initialize_performance_() throw()
{

  performance_ = new Performance (config_);

  const bool in_charm = true;
  Performance * p = performance_;
  p->new_region(perf_unknown,            "unknown");
  p->new_region(perf_simulation,         "simulation");
  p->new_region(perf_cycle,              "cycle");
  p->new_region(perf_initial,            "initial");
  p->new_region(perf_adapt_apply,        "adapt_apply");
  p->new_region(perf_adapt_apply_sync,   "adapt_apply_sync",in_charm);
  p->new_region(perf_adapt_notify,       "adapt_notify");
  p->new_region(perf_adapt_notify_sync,  "adapt_notify_sync",in_charm);
  p->new_region(perf_adapt_update,       "adapt_update");
  p->new_region(perf_adapt_update_sync,  "adapt_update_sync",in_charm);
  p->new_region(perf_adapt_end,          "adapt_end");
  p->new_region(perf_adapt_end_sync,     "adapt_end_sync",in_charm);
  p->new_region(perf_refresh_store,      "refresh_store");
  p->new_region(perf_refresh_child,      "refresh_child");
  p->new_region(perf_refresh_exit,       "refresh_exit");
  p->new_region(perf_refresh_store_sync, "refresh_store_sync",in_charm);
  p->new_region(perf_refresh_child_sync, "refresh_child_sync",in_charm);
  p->new_region(perf_refresh_exit_sync,  "refresh_exit_sync",in_charm);
  p->new_region(perf_compute,            "compute");
  p->new_region(perf_control,            "control");
  p->new_region(perf_output,             "output");
  p->new_region(perf_stopping,           "stopping");
  p->new_region(perf_block,              "block");
  p->new_region(perf_exit,               "exit");
#ifdef CONFIG_USE_GRACKLE
  p->new_region(perf_grackle,            "grackle");
#endif

  timer_.start();

#ifdef CONFIG_USE_PAPI  
  for (size_t i=0; i<config_->performance_papi_counters.size(); i++) {
    p->new_counter(counter_type_papi, 
		   config_->performance_papi_counters[i]);
  }
#endif  

#ifdef CONFIG_USE_PROJECTIONS
  int index_on = config_->performance_on_schedule_index;
  int index_off = config_->performance_off_schedule_index;
  projections_tracing_ = config_->performance_projections_on_at_start;
  if (projections_tracing_ == false) {
    traceEnd();
  }    
  if (index_on >= 0) {
    projections_schedule_on_ = Schedule::create
      ( config_->schedule_var[index_on],
	config_->schedule_type[index_on],
	config_->schedule_start[index_on],
	config_->schedule_stop[index_on],
	config_->schedule_step[index_on],
	config_->schedule_list[index_on]);
  }
  if (index_off >= 0) {
    projections_schedule_off_ = Schedule::create
      ( config_->schedule_var[index_off],
	config_->schedule_type[index_off],
	config_->schedule_start[index_off],
	config_->schedule_stop[index_off],
	config_->schedule_step[index_off],
	config_->schedule_list[index_off]);
    
  }
#endif

  p->begin();

  p->start_region(perf_simulation);

}

//----------------------------------------------------------------------

void Simulation::initialize_config_() throw()
{
  TRACE("BEGIN Simulation::initialize_config_");
  TRACE("END   Simulation::initialize_config_");
}

//----------------------------------------------------------------------

void Simulation::initialize_monitor_() throw()
{
  bool debug = config_->monitor_debug;
  int debug_mode = debug ? monitor_mode_all : monitor_mode_none;
  monitor_->set_mode("DEBUG",debug_mode);
  monitor_->set_verbose(config_->monitor_verbose);
}

//----------------------------------------------------------------------

void Simulation::initialize_data_descr_() throw()
{
  scalar_descr_long_double_ = new ScalarDescr;
  scalar_descr_double_      = new ScalarDescr;
  scalar_descr_int_         = new ScalarDescr;
  scalar_descr_long_long_   = new ScalarDescr;
  scalar_descr_sync_        = new ScalarDescr;
  scalar_descr_void_        = new ScalarDescr;
  scalar_descr_index_       = new ScalarDescr;

  //--------------------------------------------------
  // parameter: Field : list
  //--------------------------------------------------

  field_descr_ = new FieldDescr;

  // Add data fields

  for (size_t i=0; i<config_->field_list.size(); i++) {
    field_descr_->insert_permanent (config_->field_list[i]);
  }

  // Define default ghost zone depth for all fields, default value of 1

  int gx = config_->field_ghost_depth[0];
  int gy = config_->field_ghost_depth[1];
  int gz = config_->field_ghost_depth[2];

  field_descr_->set_default_ghost_depth (gx,gy,gz);

  // Default precision

  for (int i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_precision(i,config_->field_precision);
  }

  //--------------------------------------------------
  // parameter: Field : alignment
  //--------------------------------------------------

  int alignment = config_->field_alignment;

  ASSERT1 ("Simulation::initialize_data_descr_",
	  "Illegal Field:alignment parameter value %d",
	   alignment,
	   1 <= alignment );
	  
  field_descr_->set_alignment (alignment);
  
  field_descr_->set_padding (config_->field_padding);

  field_descr_->set_history (config_->field_history);

  for (int i=0; i<field_descr_->field_count(); i++) {

    std::string field_name = field_descr_->field_name(i);

    const int cx = config_->field_centering[0][i];
    const int cy = config_->field_centering[1][i];
    const int cz = config_->field_centering[2][i];

    field_descr_->set_centering(i,cx,cy,cz);

  }

  // field groups

  int num_fields = config_->field_group_list.size();
  for (int index_field=0; index_field<num_fields; index_field++) {
    std::string field = config_->field_list[index_field];
    int num_groups = config_->field_group_list[index_field].size();
   
    for (int index_group=0; index_group<num_groups; index_group++) {
      std::string group = config_->field_group_list[index_field][index_group];
      field_descr_->groups()->add(field,group);
    }
  }

  //--------------------------------------------------
  // parameter: Particle : list
  //--------------------------------------------------

  particle_descr_ = new ParticleDescr;

  // Set particle batch size
  particle_descr_->set_batch_size(config_->particle_batch_size);

  // Add particle types

  // ... first map attribute scalar type name to type_enum int
  std::map<std::string,int> type_val;
  for (int i=0; i<NUM_TYPES; i++) {
    type_val[cello::type_name[i]] = i;
  }
#ifdef CONFIG_PRECISION_SINGLE	
  type_val["default"] = type_float;
#endif  
#ifdef CONFIG_PRECISION_DOUBLE
  type_val["default"] = type_double;
#endif  
  for (size_t it=0; it<config_->particle_list.size(); it++) {

    particle_descr_->new_type (config_->particle_list[it]);

    // Add particle constants
    int nc = config_->particle_constant_name[it].size();
    for (int ic=0; ic<nc; ic++) {
      std::string name = config_->particle_constant_name[it][ic];
      int         type = type_val[config_->particle_constant_type[it][ic]];
      ASSERT3 ("Simulation::initialize_data_descr_()",
	       "Unknown Particle type \"%s\" constant \"%s\" "
	       "has unknown type \"%s\"",
	       config_->particle_list[it].c_str(),
	       name.c_str(),
	       config_->particle_attribute_type[it][ic].c_str(),
	       cello::type_is_valid(type));
      particle_descr_->new_constant(it,name,type);
      union {
	char * c;
	long long * ill;
	float * f4;
	double * f8;
	long double * f16;
	int8_t * i8;
	int16_t * i16;
	int32_t * i32;
	int64_t * i64;
      };
      c = particle_descr_->constant_value(it,ic);
      if (type == type_default) type = default_type;
      switch (type) {
      case type_single:     *f4 = config_->particle_constant_value[it][ic];
	break;
      case type_double:     *f8 = config_->particle_constant_value[it][ic];
	break;
      case type_quadruple:  *f16 = config_->particle_constant_value[it][ic];
	break;
      case type_int8:       *i8 = config_->particle_constant_value[it][ic];
	break;
      case type_int16:      *i16 = config_->particle_constant_value[it][ic];
	break;
      case type_int32:      *i32 = config_->particle_constant_value[it][ic];
	break;
      case type_int64:      *i64 = config_->particle_constant_value[it][ic];
	break;
      default:
	ERROR3 ("Simulation::initialize_data_descr_()",
		"Unrecognized type %d for particle constant %s in type %s",
		type,name.c_str(),config_->particle_list[it].c_str());
	break;
      }
    }

    // Add particle attributes
    int na = config_->particle_attribute_name[it].size();
    for (int ia=0; ia<na; ia++) {
      std::string name = config_->particle_attribute_name[it][ia];
      int type         = type_val[config_->particle_attribute_type[it][ia]];
      ASSERT3 ("Simulation::initialize_data_descr_()",
	       "Unknown Particle type \"%s\" attribute \"%s\" has unknown type \"%s\"",
	       config_->particle_list[it].c_str(),
	       name.c_str(),
	       config_->particle_attribute_type[it][ia].c_str(),
	       cello::type_is_valid(type));
      particle_descr_->new_attribute(it,name,type);
    }

    // position and velocity attributes
    particle_descr_->set_position 
      (it,
       config_->particle_attribute_position[0][it],
       config_->particle_attribute_position[1][it],
       config_->particle_attribute_position[2][it]);
    particle_descr_->set_velocity 
      (it,
       config_->particle_attribute_velocity[0][it],
       config_->particle_attribute_velocity[1][it],
       config_->particle_attribute_velocity[2][it]);
  }

  // particle groups

  int num_particles = config_->particle_group_list.size();
  for (int index_particle=0; index_particle<num_particles; index_particle++) {
    std::string particle = config_->particle_list[index_particle];
    int num_groups = config_->particle_group_list[index_particle].size();
   
    for (int index_group=0; index_group<num_groups; index_group++) {
      std::string group = config_->particle_group_list[index_particle][index_group];
      particle_descr_->groups()->add(particle,group);
    }
  }

}
//----------------------------------------------------------------------

void Simulation::initialize_hierarchy_() throw()
{
#ifdef DEBUG_SIMULATION
  CkPrintf ("%d DEBUG_SIMULATION Simulation::initialize_hierarchy_()\n",CkMyPe());
  fflush(stdout);
#endif  

  ASSERT("Simulation::initialize_hierarchy_",
	 "data must be initialized before hierarchy",
	 field_descr_ != NULL);

  //----------------------------------------------------------------------
  // Create and initialize Hierarchy
  //----------------------------------------------------------------------

  const int refinement = 2;

  hierarchy_ = factory()->create_hierarchy 
    (refinement,
     config_->mesh_min_level,
     config_->mesh_max_level);

  // Domain extents

  hierarchy_->set_lower
    (config_->domain_lower[0], 
     config_->domain_lower[1], 
     config_->domain_lower[2]);
  hierarchy_->set_upper
    (config_->domain_upper[0], 
     config_->domain_upper[1], 
     config_->domain_upper[2]);

  //----------------------------------------------------------------------
  // Create and initialize root Patch in Hierarchy
  //----------------------------------------------------------------------

  //--------------------------------------------------
  // parameter: Mesh : root_size
  // parameter: Mesh : root_blocks
  //--------------------------------------------------

  hierarchy_->set_root_size(config_->mesh_root_size[0],
			    config_->mesh_root_size[1],
			    config_->mesh_root_size[2]);

  hierarchy_->set_blocking(config_->mesh_root_blocks[0],
			   config_->mesh_root_blocks[1],
			   config_->mesh_root_blocks[2]);

  bool lp3[3] = { false, false, false };
  int index_boundary = 0;
  auto & root_blocks = config_->mesh_root_blocks;
  Boundary * boundary;
  while ( (boundary = cello::boundary(index_boundary++)) ) {
    boundary->periodicity(lp3);
  }
  int p3[3] = {lp3[0] ? root_blocks[0] : 0,
               lp3[1] ? root_blocks[1] : 0,
               lp3[2] ? root_blocks[2] : 0 };
  hierarchy_->set_periodicity (p3[0],p3[1],p3[2]);
}

//----------------------------------------------------------------------

void Simulation::initialize_balance_() throw()
{
  int index = config_->balance_schedule_index;

  schedule_balance_ = (index == -1) ? NULL : Schedule::create
    ( config_->schedule_var[index],
      config_->schedule_type[index],
      config_->schedule_start[index],
      config_->schedule_stop[index],
      config_->schedule_step[index],
      config_->schedule_list[index]);

}

//----------------------------------------------------------------------

void Simulation::initialize_block_array_() throw()
{
  bool allocate_blocks = (CkMyPe() == 0);

  // Don't allocate blocks if reading data from files

  //  bool allocate_data = ! ( config_->initial_type == "file" || 
  //			   config_->initial_type == "checkpoint" );
  bool allocate_data = true;

  if (allocate_blocks) {

    // Create the root-level blocks for level = 0
    hierarchy_->create_block_array (allocate_data);

    // Create the "sub-root" blocks if mesh_min_level < 0
    if (hierarchy_->min_level() < 0) {
      hierarchy_->create_subblock_array	(allocate_data);
    }

    hierarchy_->block_array().doneInserting();
  }
}

//----------------------------------------------------------------------

void Simulation::p_set_block_array(CProxy_Block block_array)
{
  if (CkMyPe() != 0) hierarchy_->set_block_array(block_array);
}

//----------------------------------------------------------------------

void Simulation::deallocate_() throw()
{
  delete factory_;       factory_     = 0;
  delete parameters_;    parameters_  = 0;
  delete hierarchy_;     hierarchy_ = 0;
  delete field_descr_;   field_descr_ = 0;
  delete performance_;   performance_ = 0;
}

//----------------------------------------------------------------------

const Factory * Simulation::factory() const throw()
{
  TRACE("Simulation::factory()");
  if (factory_ == NULL) factory_ = new Factory;
  return factory_;
}

//======================================================================

void Simulation::update_state(int cycle, double time, double dt, double stop) 
{
  cycle_ = cycle;
  time_  = time;
  dt_    = dt;
  stop_  = stop != 0;
}

//======================================================================

void Simulation::data_insert_block(Block * block) 
{
 
#ifdef CELLO_DEBUG
  PARALLEL_PRINTF ("%d: ++sync_output_begin_ %d %lu\n",
		   CkMyPe(),sync_output_begin_.stop(),hierarchy_->num_blocks());
#endif
  if (hierarchy_) {
    hierarchy_->insert_block(block);
    hierarchy_->increment_block_count(1,block->level());
  }
  ++sync_output_begin_;
  ++sync_output_write_;
}

//----------------------------------------------------------------------

void Simulation::data_delete_block(Block * block) 
{
  if (hierarchy_) {
    hierarchy_->delete_block(block);
    hierarchy_->increment_block_count(-1,block->level());
  }
  --sync_output_begin_;
  --sync_output_write_;
}

//----------------------------------------------------------------------

void Simulation::data_insert_particles(int64_t count)
{
  if (hierarchy_) hierarchy_->increment_particle_count(count);
}

//----------------------------------------------------------------------

void Simulation::data_delete_particles(int64_t count)
{
  if (hierarchy_) hierarchy_->increment_particle_count(-count);
}

//----------------------------------------------------------------------

void Simulation::monitor_output()
{
  monitor()-> print("", "-------------------------------------");
  monitor()-> print("Simulation", "cycle %04d", cycle_);
  monitor()-> print("Simulation", "time-sim %15.12e",time_);
  monitor()-> print("Simulation", "dt %15.12e", dt_);
  thisProxy.p_monitor_performance();
}

//----------------------------------------------------------------------

void Simulation::monitor_performance()
{
  int nr  = performance_->num_regions();
  int nc =  performance_->num_counters();

  // 0 num_sum
  // 1 num_max
  // 2 msg_coarsen
  // 3 msg_refine
  // 4 msg_refresh
  // 5 data_msg
  // 6 field_face
  // 7 particle_data
  // 8 num-particles
  // 9+ num_solver_iters
  // NL+ num-blocks-<L>
  // 10+ num_blocks_total
  // 11+ max_proc_blocks
  // 12+ max_proc_particles
  // 13+ max_node_blocks
  // 14+ max_node_particles
  // 15+ max_solver_iters
  
  const int num_solver = problem()->num_solvers();

  int n = 14 + 2*num_solver + ( hierarchy_->max_level() - hierarchy_->min_level() + 1) + nr*nc;

  
  long long * counters_region = new long long [nc];
  long long * counters_reduce = new long long [n];

  const int in = cello::index_static();
  
  int m=0;
  const int num_max = 4 + num_solver;
  counters_reduce[m++] = n - num_max - 2;
  counters_reduce[m++] = num_max;
  
  // accumulated metrics
  
  counters_reduce[m++] = MsgCoarsen::counter[in];     // 2
  counters_reduce[m++] = MsgRefine::counter[in];      // 3
  counters_reduce[m++] = MsgRefresh::counter[in];     // 4
  counters_reduce[m++] = DataMsg::counter[in];        // 5
  counters_reduce[m++] = FieldFace::counter[in];      // 6
  counters_reduce[m++] = ParticleData::counter[in];   // 7
  counters_reduce[m++] = hierarchy_->num_particles(); // 8
  for (int i=0; i<num_solver; i++) {
    counters_reduce[m++] = cello::simulation()->get_solver_num_iter(i); // 9
  }

  const int min_level = hierarchy_->min_level();

  int num_blocks_total = 0;
  for (int i=min_level; i<=hierarchy_->max_level(); i++) {
    num_blocks_total +=  hierarchy_->num_blocks(i);
    counters_reduce[m++] = hierarchy_->num_blocks(i); // NL
  }
  counters_reduce[m++] = num_blocks_total;            // 10  num_blocks_total
  
  // performance region counters
  for (int ir = 0; ir < nr; ir++) {
    performance_->region_counters(ir,counters_region);
    for (int ic = 0; ic < nc; ic++) {
      counters_reduce[m++] = counters_region[ic];
    }
  }

  // maximum metrics
  
  counters_reduce[m++] = num_blocks_total;            // 11  max_proc_blocks
  counters_reduce[m++] = hierarchy_->num_particles(); // 12  max_proc_particles
  counters_reduce[m++] = Hierarchy::num_blocks_node;  // 13  max_node_blocks
  counters_reduce[m++] = Hierarchy::num_particles_node;// 14 max_node_particles
  for (int i=0; i<num_solver; i++) {
    counters_reduce[m++] = cello::simulation()->get_solver_max_iter(i); // 15 max_node_particles
  }

  ASSERT2("Simulation::monitor_performance()",
	  "Actual array length %d != expected array length %d", m,n,
	  (m == n) );
  
  // --------------------------------------------------
#ifdef TRACE_CONTRIBUTE
  CkPrintf ("%s:%d DEBUG_CONTRIBUTE\n",__FILE__,__LINE__); fflush(stdout);
#endif  

  contribute
    (n*sizeof(long long),
     counters_reduce,
     r_reduce_performance_type,
     CkCallback (CkIndex_Simulation::r_monitor_performance_reduce(NULL), 
		 thisProxy));
  // --------------------------------------------------

  delete [] counters_reduce;
  delete [] counters_region;

}

//----------------------------------------------------------------------

void Simulation::r_monitor_performance_reduce(CkReductionMsg * msg)
{
  long long * counters_reduce = (long long *)msg->getData();

  int index_region_cycle = performance_->region_index("cycle");

  int m = 0;
  const int num_sum = counters_reduce[m++];             // 0
  const int num_max = counters_reduce[m++];             // 1
  const long long msg_coarsen = counters_reduce[m++];   // 2
  const long long msg_refine  = counters_reduce[m++];   // 3
  const long long msg_refresh = counters_reduce[m++];   // 4
  const long long data_msg    = counters_reduce[m++];   // 5
  const long long field_face  = counters_reduce[m++];   // 6
  const long long particle_data = counters_reduce[m++]; // 7
  const long long num_particles = counters_reduce[m++]; // 8

  const int num_solver = problem()->num_solvers();
  for (int i=0; i<num_solver; i++) {
    const long long num_solver_iter = counters_reduce[m++]; // 15
    monitor()->print ("Performance","solver num-%s-iter %lld",
                      problem()->solver(i)->name().c_str(),
                      num_solver_iter);
  }

  monitor()->print("Performance","counter num-msg-coarsen %lld", msg_coarsen);
  monitor()->print("Performance","counter num-msg-refine %lld", msg_refine);
  monitor()->print("Performance","counter num-msg-refresh %lld", msg_refresh);
  monitor()->print("Performance","counter num-data-msg %lld", data_msg);
  monitor()->print("Performance","counter num-field-face %lld", field_face);
  monitor()->print("Performance","counter num-particle-data %lld", particle_data);

  monitor()->print("Performance","simulation num-particles total %lld",
		   num_particles);

  // compute total blocks and leaf blocks
  long long num_total_blocks = 0;
  long long num_leaf_blocks = 0;
  for (int i=hierarchy_->min_level(); i<=hierarchy_->max_level(); i++) {
    const long long num_blocks_level = counters_reduce[m++]; // NL
    monitor()->print("performance","simulation num-blocks-level %d %lld",
		     i,num_blocks_level);
    num_total_blocks += num_blocks_level;
    // compute leaf blocks given number of blocks per level
    // (NOTE: num_blocks_level (i>0) is evenly divisible by num_children
    if (i==0) {
      num_leaf_blocks = num_blocks_level;
    } else if (i>0) {
      num_leaf_blocks +=
	(num_blocks_level - num_blocks_level/cello::num_children());
    }
  }

  monitor()->print
    ("Performance","simulation num-leaf-blocks %lld",  num_leaf_blocks);
  monitor()->print
    ("Performance","simulation num-total-blocks %lld", num_total_blocks);

  const long long num_blocks_total   = counters_reduce[m++]; // 10

  if (num_total_blocks != num_blocks_total) {
    WARNING2 ("Simulation::r_monitor_performance_reduce()",
	      "num_blocks_total %lld does not match computed value %lld",
	      num_total_blocks,num_blocks_total);
  }
  
  const int num_regions  = performance_->num_regions();
  const int num_counters =  performance_->num_counters();

  for (int ir = 0; ir < num_regions; ir++) {
    for (int ic = 0; ic < num_counters; ic++, m++) {
      bool do_print =
	(ir != perf_unknown) && (
	(performance_->counter_type(ic) != counter_type_abs) ||
	(ir == index_region_cycle));
      if (do_print) {
	monitor()->print("Performance","%s %s %lld",
			performance_->region_name(ir).c_str(),
			performance_->counter_name(ic).c_str(),
			 counters_reduce[m]);
      }
      
    }
  }

  const long long max_proc_blocks    = counters_reduce[m++]; // 11
  const long long max_proc_particles = counters_reduce[m++]; // 12
  const long long max_node_blocks    = counters_reduce[m++]; // 13
  const long long max_node_particles = counters_reduce[m++]; // 14

  for (int i=0; i<num_solver; i++) {
    const long long max_solver_iters       = counters_reduce[m++]; // 15
    monitor()->print ("Performance","solver max-%s-iter %lld",
                      problem()->solver(i)->name().c_str(),
                      max_solver_iters);
  }
  cello::simulation()->clear_solver_iter(); // clear it for the next solve

  
  monitor()->print
    ("Performance","simulation max-proc-blocks %lld",  max_proc_blocks);
  monitor()->print
    ("Performance","simulation max-node-blocks %lld",  max_node_blocks);
  monitor()->print
    ("Performance","simulation max-proc-particles %lld", max_proc_particles);
  monitor()->print
    ("Performance","simulation max-node-particles %lld", max_node_particles);

  const double avg_proc_blocks = 1.0*num_blocks_total/CkNumPes();
  const double avg_node_blocks = 1.0*num_blocks_total/CkNumNodes();


  // monitor()->print
  //   ("Performance","simulation balance-blocks-core %f",
  //    100.0*(max_proc_blocks / avg_proc_blocks - 1.0 ));
  // monitor()->print
  //   ("Performance","simulation balance-blocks-node %f",
  //    100.0*(max_node_blocks / avg_node_blocks - 1.0 ));

  monitor()->print
    ("Performance","simulation balance-eff-blocks-core %f (%.0f/%lld)",
     avg_proc_blocks / max_proc_blocks,
     avg_proc_blocks, max_proc_blocks);
  monitor()->print
    ("Performance","simulation balance-eff-blocks-node %f (%.0f/%lld)",
     avg_node_blocks / max_node_blocks,
     avg_node_blocks, max_node_blocks);
  

  if (num_particles > 0) {
    const double avg_proc_particles = 1.0*num_particles/CkNumPes();
    const double avg_node_particles = 1.0*num_particles/CkNumNodes();
    monitor()->print
      ("Performance","simulation balance-eff-particles-core %f (%.0f/%lld)",
       avg_proc_particles / max_proc_particles,
       avg_proc_particles , max_proc_particles );
    monitor()->print
      ("Performance","simulation balance-eff-particles-node %f (%.0f/%lld)",
       avg_node_particles / max_node_particles,
       avg_node_particles , max_node_particles );
  }
  
  ASSERT3("Simulation::monitor_performance()",
	  "Actual array length %d != expected array length 2 + %d + %d",
	  m,num_sum,num_max,
	  (m == 2+num_sum+num_max) );

  delete msg;

  Memory::instance()->reset_high();

}
