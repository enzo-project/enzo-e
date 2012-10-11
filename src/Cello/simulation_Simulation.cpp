// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "main.hpp"

#include "simulation.hpp"
#include "charm_simulation.hpp"

bool config_performance = false;

Simulation::Simulation
(
 const char *   parameter_file,
#ifdef CONFIG_USE_CHARM
 int            n,
#endif
 const GroupProcess * group_process
 ) throw()
/// Initialize the Simulation object
: factory_(0),
  parameters_(0),
  parameter_file_(parameter_file),
  group_process_((GroupProcess *)group_process),
  is_group_process_new_(false),
#ifdef CONFIG_USE_CHARM
  patch_loop_(0),
#endif
  dimension_(0),
  cycle_(0),
  level_(0),
  time_(0.0),
  dt_(0),
  stop_(false),
  performance_simulation_(0),
  performance_cycle_(0),
  performance_curr_(0),
  lcaperf_(0),
  monitor_(0),
  hierarchy_(0),
  field_descr_(0)
{

  if (!group_process_) {
    group_process_ = GroupProcess::create();
    is_group_process_new_ = true;
  }

  monitor_ = Monitor::instance();
  monitor_->set_process_rank(group_process_->rank());
  monitor_->set_active(group_process_->is_root());


  num_perf_ = 1;

#ifdef CONFIG_USE_PAPI
  num_perf_ += 4;
#endif

  perf_val_ = new double[num_perf_];
  perf_min_ = new double[num_perf_];
  perf_max_ = new double[num_perf_];
  perf_sum_ = new double[num_perf_];

  performance_simulation_ = new Performance;
  performance_cycle_      = new Performance;


  lcaperf_ = new LcaPerf (group_process_->rank(), group_process_->size());

  lcaperf_->initialize();

  lcaperf_->new_region("simulation");
  lcaperf_->new_attribute("cycle",LCAP_INT);
  lcaperf_->new_attribute("level",LCAP_INT);

  lcaperf_->begin();

  parameters_ = new Parameters(parameter_file,monitor_);

  config_.read(parameters_);
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

Simulation::Simulation()
  : patch_loop_(0),
    lcaperf_(0)
{ TRACE("Simulation()"); }

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
void Simulation::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  CBase_Simulation::pup(p);
  bool up = p.isUnpacking();

  p | factory_; // PUP::able

  // if (up) parameters_ = new Parameters;
  // p | * parameters_;

  p | config_;

  p | parameter_file_;

  if (up) group_process_ = GroupProcess::create();

  p | is_group_process_new_;
  p | patch_loop_;
  p | dimension_; 
  p | cycle_;
  p | time_;
  p | dt_;
  p | stop_;

  if (up) problem_ = new Problem;
  p | * problem_;

  if (up) performance_simulation_ = new Performance;
  p | * performance_simulation_;

  if (up) performance_cycle_ = new Performance;
  p | * performance_cycle_;

  if (up) performance_curr_ = new Performance;
  p | * performance_curr_;

  if (up) lcaperf_ = new LcaPerf;
  p | *lcaperf_;

  p | num_perf_;
  if (up) perf_val_ = new double [num_perf_];
  PUParray(p,perf_val_,num_perf_);
  if (up) perf_min_ = new double [num_perf_];
  PUParray(p,perf_min_,num_perf_);
  if (up) perf_max_ = new double [num_perf_];
  PUParray(p,perf_max_,num_perf_);
  if (up) perf_sum_ = new double [num_perf_];
  PUParray(p,perf_sum_,num_perf_);

  if (up) monitor_ = Monitor::instance();

  if (up) hierarchy_ = new Hierarchy;
  p | *hierarchy_;

  if (up) field_descr_ = new FieldDescr;
  p | *field_descr_;
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

Simulation::Simulation (CkMigrateMessage *m)
  : CBase_Simulation(m),
    patch_loop_(0),
    lcaperf_(0)
{ TRACE("Simulation(Ck`MigrateMessage)"); }

#endif

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Simulation::initialize() throw()
{
  initialize_monitor_();

  initialize_simulation_();


  initialize_data_descr_();

  problem_->initialize_boundary(&config_);
  problem_->initialize_initial (&config_,parameters_,group_process_);
  problem_->initialize_stopping(&config_);
  problem_->initialize_timestep(&config_);
  problem_->initialize_output  (&config_,field_descr_,group_process_,factory());
  problem_->initialize_method  (&config_);

  initialize_hierarchy_();

}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  DEBUG0;

  lcaperf_->stop("simulation");
  lcaperf_->end();
  lcaperf_->finalize();

  performance_simulation_->stop();
  performance_cycle_->stop();
}

//======================================================================

void Simulation::initialize_simulation_() throw()
{

  dimension_ = config_.mesh_root_rank;
  
  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be specified",
	  dimension_ != 0);
  
  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be 1, 2, or 3",
	  (1 <= dimension_) && (dimension_ <= 3));

  cycle_ = config_.initial_cycle;
  time_  = config_.initial_time;
  dt_ = 0;

  // Initialize Performance

  lcaperf_->attribute("cycle",&cycle_,LCAP_INT);
  lcaperf_->attribute("level",&level_,LCAP_INT);
  lcaperf_->start("simulation");

  performance_simulation_->start();
  performance_cycle_->start();

}

//----------------------------------------------------------------------

void Simulation::initialize_monitor_() throw()
{
  bool debug = config_.monitor_debug;
  monitor_->set_active("DEBUG",debug);
 
}

//----------------------------------------------------------------------

void Simulation::initialize_data_descr_() throw()
{

  field_descr_ = new FieldDescr;

  //--------------------------------------------------
  // parameter: Field : fields
  //--------------------------------------------------

  // Add data fields

  for (size_t i=0; i<config_.field_fields.size(); i++) {
    field_descr_->insert_field (config_.field_fields[i]);
  }

  // Define default ghost zone depth for all fields, default value of 1

  int gx = config_.field_ghosts[0];
  int gy = config_.field_ghosts[1];
  int gz = config_.field_ghosts[2];

  for (int i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_ghosts (i,gx,gy,gz);
  }

  // Set face dimensions to refresh

  field_descr_->set_refresh_face(2,config_.field_refresh_faces);
  field_descr_->set_refresh_face(1,config_.field_refresh_edges);
  field_descr_->set_refresh_face(0,config_.field_refresh_corners);
  
  // Default precision

  for (int i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_precision(i,config_.field_precision);
  }

  //--------------------------------------------------
  // parameter: Field : alignment
  //--------------------------------------------------

  int alignment = config_.field_alignment;

  ASSERT1 ("Simulation::initialize_data_descr_",
	  "Illegal Field:alignment parameter value %d",
	   alignment,
	   1 <= alignment );
	  
  field_descr_->set_alignment (alignment);
  
  field_descr_->set_padding (config_.field_padding);


  for (int i=0; i<field_descr_->field_count(); i++) {

    std::string field_name = field_descr_->field_name(i);

    bool cx = config_.field_centering[0][i];
    bool cy = config_.field_centering[1][i];
    bool cz = config_.field_centering[2][i];

    field_descr_->set_centering(i,cx,cy,cz);

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
  hierarchy_ = factory()->create_hierarchy (dimension_,refinement);

  // Domain extents

  hierarchy_->set_lower
    (config_.domain_lower[0], 
     config_.domain_lower[1], 
     config_.domain_lower[2]);
  hierarchy_->set_upper
    (config_.domain_upper[0], 
     config_.domain_upper[1], 
     config_.domain_upper[2]);

  //----------------------------------------------------------------------
  // Create and initialize root Patch in Hierarchy
  //----------------------------------------------------------------------

  //--------------------------------------------------
  // parameter: Mesh : root_size
  // parameter: Mesh : root_blocks
  //--------------------------------------------------

  hierarchy_->set_root_size(config_.mesh_root_size[0],
			    config_.mesh_root_size[1],
			    config_.mesh_root_size[2]);

  // Don't allocate blocks if reading data from files

  bool allocate_blocks = ! ( config_.initial_type == "file" || 
			     config_.initial_type == "restart" );

#ifdef CONFIG_USE_CHARM
  // Distributed patches in Charm: only allocate on root processor
  ++patch_loop_.stop();
  if (group_process()->is_root())
#endif
    {
      hierarchy_->create_root_patch
	(field_descr_,
	 config_.mesh_root_size[0],
	 config_.mesh_root_size[1],
	 config_.mesh_root_size[2],
	 config_.mesh_root_blocks[0],
	 config_.mesh_root_blocks[1],
	 config_.mesh_root_blocks[2],
	 allocate_blocks);
    }

}

//----------------------------------------------------------------------

void Simulation::deallocate_() throw()
{
  delete factory_;       factory_     = 0;
  delete parameters_;    parameters_  = 0;
  delete performance_simulation_; performance_simulation_ = 0;
  delete performance_cycle_;      performance_cycle_ = 0;
  delete lcaperf_; lcaperf_ = 0;
  delete [] perf_val_;
  delete [] perf_min_;
  delete [] perf_max_;
  delete [] perf_sum_;
  if (is_group_process_new_)
    { delete group_process_; group_process_ = 0; }
  delete hierarchy_;     hierarchy_ = 0;
  delete field_descr_;   field_descr_ = 0;
}

//----------------------------------------------------------------------

const Factory * Simulation::factory() const throw()
{
  DEBUG("Simulation::factory()");
  if (factory_ == NULL) factory_ = new Factory;
  return factory_;
}

//======================================================================


#ifdef CONFIG_USE_CHARM

#endif

#ifndef CONFIG_USE_CHARM

//----------------------------------------------------------------------
// NOT CHARM
//----------------------------------------------------------------------


#endif

//----------------------------------------------------------------------

void Simulation::update_state(int cycle, double time, double dt, double stop) 
{
  DEBUG4 ("Simulation::update_state cycle %d time %f dt %f stop %f",
	  cycle,time,dt,stop);
 
  cycle_ = cycle;
  time_  = time;
  dt_    = dt;
  stop_  = stop;
}

//----------------------------------------------------------------------

void Simulation::monitor_output()
{

  monitor_->  print("", "-------------------------------------");

  monitor_-> print("Simulation", "cycle %04d", cycle_);
  monitor_-> print("Simulation", "time-sim %15.12f",time_);
  monitor_-> print("Simulation", "dt %15.12g", dt_);

  Memory * memory = Memory::instance();

  if (memory->is_active()) {
    monitor_->print("Memory","bytes-curr %lld", memory->bytes());
    monitor_->print("Memory","bytes-high %lld", memory->bytes_high());

    memory->reset_high();
  }

  if (config_performance) {
    performance_output(performance_cycle_);
  } else {
#ifdef CONFIG_USE_CHARM
    ((SimulationCharm *) this)->c_compute();
#endif
  }

}


//----------------------------------------------------------------------

void Simulation::performance_output(Performance * performance)
{
  lcaperf_->stop("simulation");
  lcaperf_->print();
  lcaperf_->attribute("cycle",&cycle_,LCAP_INT);
  lcaperf_->start("simulation");

  performance_curr_ = performance;

  size_t i = 0;

  // Real time
  perf_val_[i++] = performance->time();

#ifdef CONFIG_USE_PAPI

  Papi * papi = performance->papi();

  papi->update();

  // PAPI real time
  perf_val_[i++]= papi->time_real();

  // PAPI proc time
  double time_real = perf_val_[i++]= papi->time_proc();

  // PAPI gflop count
  double gflop_count = perf_val_[i++]= papi->flop_count()*1e-9;

  // PAPI gflop rate
  perf_val_[i++]= gflop_count / time_real;

#endif

#ifdef CONFIG_USE_CHARM

  // Save the performance object

  // First reduce minimum values

  CkCallback callback (CkIndex_SimulationCharm::p_performance_min(NULL),
		       thisProxy);

  contribute( num_perf_*sizeof(double), perf_val_, 
	      CkReduction::min_double, callback);
#else

  Reduce * reduce = group_process()->create_reduce();

  for (size_t i = 0; i < num_perf_; i++) {
    perf_min_[i] = reduce->reduce_double(perf_val_[i],reduce_op_min);
    perf_max_[i] = reduce->reduce_double(perf_val_[i],reduce_op_max);
    perf_sum_[i] = reduce->reduce_double(perf_val_[i],reduce_op_sum);
  }

  delete reduce; reduce = 0;

  output_performance_();

#endif

}

//----------------------------------------------------------------------

void Simulation::output_performance_()
{
  DEBUG("Simulation::output_performance");
  int i = 0;
  int np = group_process()->size();

  std::string region;

  if (performance_curr_ == performance_simulation_) {
    region = "total";
  } else if (performance_curr_ == performance_cycle_) {
    region = "cycle";
  } else {
    ERROR1 ("Simulation::output_performance_",
	    "Illegal performance_curr_ pointer %p",
	    performance_curr_);
  }

  monitor_->print ("Performance","%s time-real        %f %f %f",
		   region.c_str(),perf_min_[i],perf_sum_[i]/np,perf_max_[i]);
  ++i;

#ifdef CONFIG_USE_PAPI

  monitor_->print ("Performance","%s time-real-papi   %f %f %f", 
		   region.c_str(), perf_min_[i],perf_sum_[i]/np,perf_max_[i]);
  ++i;
  monitor_->print ("Performance","%s time-proc-papi   %f %f %f",
		   region.c_str(), perf_min_[i],perf_sum_[i]/np,perf_max_[i]);
  ++i;
  monitor_->print ("Performance","%s gflop-count-papi %f %f %f",
		   region.c_str(), perf_min_[i],perf_sum_[i]/np,perf_max_[i]);
  ++i;
  monitor_->print ("Performance","%s gflop-rate-papi  %f %f %f",
		   region.c_str(), perf_min_[i],perf_sum_[i]/np,perf_max_[i]);
#endif

#ifdef CONFIG_USE_CHARM
  if (performance_curr_ == performance_cycle_) {
    ((SimulationCharm *) this)->c_compute();
  // } else {
  //   DEBUG("Calling p_exit");
  //   proxy_main.p_exit(CkNumPes());
  }
#endif

  // monitor_->set_active(save_active);

}
//======================================================================
