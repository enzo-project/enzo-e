// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "main.hpp"

#include "simulation.hpp"
#include "charm_simulation.hpp"

Simulation::Simulation
(
 const char *   parameter_file,
 int            n
 )
/// Initialize the Simulation object
: factory_(0),
  parameters_(0),
  parameter_file_(parameter_file),
  rank_(0),
  cycle_(0),
  time_(0.0),
  dt_(0),
  stop_(false),
  phase_(phase_unknown),
  config_(0),
  problem_(0),
  timer_(),
  performance_(NULL),
  performance_name_(""),
  performance_stride_(1),
  // projections_tracing_(1),
  schedule_balance_(0),
  monitor_(0),
  hierarchy_(0),
  field_descr_(0),
  particle_descr_(0)
{
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
{ TRACE("Simulation()"); }

//----------------------------------------------------------------------

void Simulation::pup (PUP::er &p)
{
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
  p | time_;
  p | dt_;
  p | stop_;
  p | phase_;

  if (up) problem_ = new Problem;
  p | * problem_;

  if (up) performance_ = new Performance;
  p | *performance_;

  p | performance_name_;
  p | performance_stride_;

  // p | projections_tracing_;
  // if (up) projections_schedule_on_ = new Schedule;
  // p | *projections_schedule_on_;
  // if (up) projections_schedule_off_ = new Schedule;
  // p | *projections_schedule_off_;

  if (up) monitor_ = Monitor::instance();
  p | *monitor_;

  if (up) hierarchy_ = new Hierarchy;
  p | *hierarchy_;

  if (up) field_descr_ = new FieldDescr;
  p | *field_descr_;

  if (up) particle_descr_ = new ParticleDescr;
  p | *particle_descr_;

  if (up && (phase_ == phase_restart)) {
    monitor_->print ("Simulation","restarting");
  }

  p | sync_output_begin_;
  p | sync_output_write_;

  if (up) sync_output_begin_.set_stop(0);
  if (up) sync_output_write_.set_stop(0);

  p | schedule_balance_;
}

//----------------------------------------------------------------------

Simulation::Simulation (CkMigrateMessage *m)
  : CBase_Simulation(m)
{ TRACE("Simulation(CkMigrateMessage)"); }

//----------------------------------------------------------------------

Simulation::~Simulation()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  TRACE0;

  performance_->stop_region(perf_simulation);

  performance_->end();


}

//======================================================================

void Simulation::initialize_simulation_() throw()
{

  rank_ = config_->mesh_root_rank;
  
  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be specified",
	  rank_ != 0);
  
  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be 1, 2, or 3",
	  (1 <= rank_) && (rank_ <= 3));

  cycle_ = config_->initial_cycle;
  time_  = config_->initial_time;
  dt_ = 0;
}

//----------------------------------------------------------------------

void Simulation::initialize_memory_() throw()
{
  Memory * memory = Memory::instance();
  if (memory) memory->set_active(config_->memory_active);
  
}
//----------------------------------------------------------------------

void Simulation::initialize_performance_() throw()
{

  performance_ = new Performance (config_);

  performance_->new_region(perf_unknown,    "unknown");
  performance_->new_region(perf_simulation, "simulation");
  performance_->new_region(perf_cycle,      "cycle");
  performance_->new_region(perf_initial,    "initial");
  performance_->new_region(perf_adapt,      "adapt");
  performance_->new_region(perf_refresh,    "refresh");
  performance_->new_region(perf_compute,    "compute");
  performance_->new_region(perf_output,     "output");
  performance_->new_region(perf_stopping,   "stopping");

  performance_name_   = config_->performance_name;
  performance_stride_ = config_->performance_stride;

  timer_.start();

  for (size_t i=0; i<config_->performance_papi_counters.size(); i++) {
    performance_->new_counter(counter_type_papi, 
			      config_->performance_papi_counters[i]);
  }

  performance_->begin();

  performance_->start_region(perf_simulation);

}

//----------------------------------------------------------------------

void Simulation::initialize_config_() throw()
{
  TRACE("BEGIN Simulation::initialize_config_");
  if (config_ == NULL) {
    config_ = new Config;
    TRACE("Simulation::initialize_config_ calling Config::read()");
    config_->read(parameters_);
  }
  if (CkMyPe() == 0) {
    parameters_->write("parameters.out");
  }
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

  for (int i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_ghost_depth (i,gx,gy,gz);
  }

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

  for (size_t it=0; it<config_->particle_list.size(); it++) {

    particle_descr_->new_type (config_->particle_list[it]);

    // Add particle constants
    int nc = config_->particle_constant_name[it].size();
    for (int ic=0; ic<nc; ic++) {
      std::string name = config_->particle_constant_name[it][ic];
      int         type = type_val[config_->particle_constant_type[it][ic]];
      particle_descr_->new_constant(it,name,type);
      union {
	char * c;
	long long * ill;
	float * f;
	double * d;
	long double * ld;
	int8_t * i8;
	int16_t * i16;
	int32_t * i32;
	int64_t * i64;
      };
      c = particle_descr_->constant_value(it,ic);
      if (type == type_default) type = default_type;
      switch (type) {
      case type_single:     *f = config_->particle_constant_value[it][ic];
	break;
      case type_double:     *d = config_->particle_constant_value[it][ic];
	break;
      case type_quadruple:  *ld = config_->particle_constant_value[it][ic];
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

  ASSERT("Simulation::initialize_hierarchy_",
	 "data must be initialized before hierarchy",
	 field_descr_ != NULL);

  //----------------------------------------------------------------------
  // Create and initialize Hierarchy
  //----------------------------------------------------------------------

  const int refinement = 2;
  hierarchy_ = factory()->create_hierarchy 
    (rank_,refinement,config_->mesh_max_level);

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

#ifdef TEMP_BALANCE_MANUAL
  if (schedule_balance_) {
    TurnManualLBOn();
  }
#endif

}

//----------------------------------------------------------------------

void Simulation::initialize_forest_() throw()
{
  bool allocate_blocks = (CkMyPe() == 0);

  // Don't allocate blocks if reading data from files

  //  bool allocate_data = ! ( config_->initial_type == "file" || 
  //			   config_->initial_type == "checkpoint" );
  bool allocate_data = true;

  if (allocate_blocks) {

    // Create the root-level blocks for level = 0
    hierarchy_->create_forest
      (field_descr_,
       allocate_data);

    // Create the "sub-root" blocks if mesh_min_level < 0
    if (config_->mesh_min_level < 0) {
      hierarchy_->create_subforest
	(field_descr_,
	 allocate_data,
	 config_->mesh_min_level);
    }
    hierarchy_->block_array()->doneInserting();


  }
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

void Simulation::monitor_insert_block(int count) 
{
 
#ifdef CELLO_DEBUG
  PARALLEL_PRINTF ("%d: ++sync_output_begin_ %d %d\n",
		   CkMyPe(),sync_output_begin_.stop(),hierarchy()->num_blocks());
#endif
  hierarchy()->increment_block_count(count);
  sync_output_begin_ += count;
  sync_output_write_ += count;
}

//----------------------------------------------------------------------

void Simulation::monitor_delete_block(int count) 
{
  hierarchy()->increment_block_count(-count);
  sync_output_begin_ -= count;
  sync_output_write_ -= count;
}

//----------------------------------------------------------------------

void Simulation::monitor_insert_zones(int64_t count_total, int64_t count_real) 
{
 
  ASSERT2 ("Simulation::monitor_insert_zones()",
	   "Total number of zones %ld must be no larger than "
	   "number of real zones %ld",
	   count_total,count_real,
	   count_total >= count_real);
	   
  hierarchy()->increment_total_zone_count(count_total);
  hierarchy()->increment_real_zone_count (count_real);
}

//----------------------------------------------------------------------

void Simulation::monitor_delete_zones(int64_t count_total, int64_t count_real)
{
  ASSERT2 ("Simulation::monitor_insert_zones()",
	   "Total number of zones %ld must be no larger than "
	   "number of real zones %ld",
	   count_total,count_real,
	   count_total >= count_real);
	   
  hierarchy()->increment_total_zone_count(-count_total);
  hierarchy()->increment_real_zone_count(-count_real);
}

//----------------------------------------------------------------------

void Simulation::monitor_insert_particles(int64_t count)
{
  hierarchy()->increment_particle_count(count);
}

//----------------------------------------------------------------------

void Simulation::monitor_delete_particles(int64_t count)
{
  hierarchy()->increment_particle_count(-count);
}

//----------------------------------------------------------------------

void Simulation::p_monitor()
{
  monitor()-> print("", "-------------------------------------");
  monitor()-> print("Simulation", "cycle %04d", cycle_);
  monitor()-> print("Simulation", "time-sim %15.12e",time_);
  monitor()-> print("Simulation", "dt %15.12e", dt_);
  proxy_simulation.p_monitor_performance();
}

//----------------------------------------------------------------------

void Simulation::monitor_performance()
{
  int nr  = performance_->num_regions();
  int nc =  performance_->num_counters();

  int n = nr * nc + 4;

  long long * counters_long_long = new long long [nc];
  long *      counters_long = new long [n];

  for (int ir = 0; ir < nr; ir++) {
    performance_->region_counters(ir,counters_long_long);
    for (int ic = 0; ic < nc; ic++) {
      int ia = ir+nr*ic;
      counters_long[ia] = (long) counters_long_long[ic];
    }
  }

  // number of Blocks
  counters_long[n-4] = hierarchy()->num_blocks(); 
  // number of particles
  counters_long[n-3] = hierarchy()->num_particles();
  // number of real zones
  counters_long[n-2] = hierarchy()->num_zones_real();
  // number of zones total
  counters_long[n-1] = hierarchy()->num_zones_total(); 

  // --------------------------------------------------
  CkCallback callback (CkIndex_Simulation::r_monitor_performance(NULL), 
		       thisProxy);
  contribute (n*sizeof(long), counters_long,CkReduction::sum_long,callback);
  // --------------------------------------------------

  delete [] counters_long;
  delete [] counters_long_long;

}

//----------------------------------------------------------------------

void Simulation::r_monitor_performance(CkReductionMsg * msg)
{
  int nr  = performance_->num_regions();
  int nc =  performance_->num_counters();

  int n = nr * nc + 4;

  long *      counters_long = (long * )msg->getData();

  int index_region_cycle = performance_->region_index("cycle");

  for (int ir = 0; ir < nr; ir++) {
    for (int ic = 0; ic < nc; ic++) {
      int ia = ir+nr*ic;
      bool do_print = 
	(performance_->counter_type(ic) != counter_type_abs) ||
	(ir == index_region_cycle);
      if (do_print) {
	monitor()->print("Performance","%s %s %ld",
			performance_->region_name(ir).c_str(),
			performance_->counter_name(ic).c_str(),
			counters_long[ia]);
      }
    }
  }

  monitor()->print("Performance","simulation num-blocks %d",
		  counters_long[n-4]);
  monitor()->print("Performance","simulation num-particles %d",
		  counters_long[n-3]);
  monitor()->print("Performance","simulation num-zones-real %d",
		  counters_long[n-2]);
  monitor()->print("Performance","simulation num-zones-total %d",
		  counters_long[n-1]);

  Memory::instance()->reset_high();

  delete msg;

}
