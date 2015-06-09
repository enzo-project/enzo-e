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
 ) throw()
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
  monitor_(0),
  hierarchy_(0),
  field_descr_(0)
{
  debug_open();

  monitor_ = Monitor::instance();
#ifdef CELLO_DEBUG
  monitor_->set_active(true);
#else
  monitor_->set_active(CkMyPe() == 0);
#endif

  // Read in parameters

  parameters_ = new Parameters(parameter_file,monitor_);

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

  if (up && (phase_ == phase_restart)) {
    monitor_->print ("Simulation","restarting");
  }

  p | sync_output_begin_;
  p | sync_output_write_;

  if (up) sync_output_begin_.set_stop(0);
  if (up) sync_output_write_.set_stop(0);
}

//----------------------------------------------------------------------

Simulation::Simulation (CkMigrateMessage *m)
  : CBase_Simulation(m)
{ TRACE("Simulation(CkMigrateMessage)"); }

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
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
  monitor_->set_active("DEBUG",debug);
  monitor_->set_verbose(config_->monitor_verbose);
}

//----------------------------------------------------------------------

void Simulation::initialize_data_descr_() throw()
{

  field_descr_ = new FieldDescr;

  //--------------------------------------------------
  // parameter: Field : fields
  //--------------------------------------------------

  // Add data fields

  for (size_t i=0; i<config_->field_list.size(); i++) {
    field_descr_->insert_permanent (config_->field_list[i]);
  }

  // Define default ghost zone depth for all fields, default value of 1

  int gx = config_->field_ghosts[0];
  int gy = config_->field_ghosts[1];
  int gz = config_->field_ghosts[2];

  for (int i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_ghosts (i,gx,gy,gz);
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

  // groups

  int num_fields = config_->field_group_list.size();
  for (int index_field=0; index_field<num_fields; index_field++) {
    std::string field = config_->field_list[index_field];
    int num_groups = config_->field_group_list[index_field].size();
   
    for (int index_group=0; index_group<num_groups; index_group++) {
      std::string group = config_->field_group_list[index_field][index_group];
      field_descr_->groups()->add(field,group);
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
    (rank_,refinement, 0, CkNumPes());

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

void Simulation::initialize_forest_() throw()
{

  bool allocate_blocks = (CkMyPe() == 0);

  // Don't allocate blocks if reading data from files

  bool allocate_data = ! ( config_->initial_type == "file" || 
			   config_->initial_type == "checkpoint" );

  if (allocate_blocks) {
    hierarchy_->create_forest
      (field_descr_,
       allocate_data);
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

void Simulation::insert_block() 
{
 
#ifdef CELLO_DEBUG
  PARALLEL_PRINTF ("%d: ++sync_output_begin_ %d %d\n",
		   CkMyPe(),sync_output_begin_.stop(),hierarchy()->num_blocks());
#endif
  hierarchy()->increment_block_count(1);
  ++sync_output_begin_;
  ++sync_output_write_;
}

//----------------------------------------------------------------------

void Simulation::delete_block() 
{
  hierarchy()->increment_block_count(-1);
  --sync_output_begin_;
  --sync_output_write_;
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

  int n = nr * nc + 1;

  long long * counters_long_long = new long long [nc];
  long *      counters_long = new long [n];

  for (int ir = 0; ir < nr; ir++) {
    performance_->region_counters(ir,counters_long_long);
    for (int ic = 0; ic < nc; ic++) {
      int index_counter = ir+nr*ic;
      counters_long[index_counter] = 
	(long) counters_long_long[ic];
    }
  }

  counters_long[n-1] = hierarchy()->num_blocks(); // number of Blocks

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

  int n = nr * nc + 1;

  long *      counters_long = (long * )msg->getData();

  int index_region_cycle = performance_->region_index("cycle");

  for (int ir = 0; ir < nr; ir++) {
    for (int ic = 0; ic < nc; ic++) {
      int index_counter = ir+nr*ic;
      bool do_print = 
	(performance_->counter_type(ic) != counter_type_abs) ||
	(ir == index_region_cycle);
      if (do_print) {
	monitor()->print("Performance","%s %s %ld",
			performance_->region_name(ir).c_str(),
			performance_->counter_name(ic).c_str(),
			counters_long[index_counter]);
      }
    }
  }

  monitor()->print("Performance","simulation num-blocks %d",
		  counters_long[n-1]);

  Memory::instance()->reset_high();

  delete msg;

}

