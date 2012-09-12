// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "main.hpp"

#include "simulation.hpp"
#include "charm_simulation.hpp"

bool config_performance = true;

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


  //  lcaperf_ = new lca::LcaPerf;

  lcaperf_.initialize();
  lcaperf_.begin();
  lcaperf_.new_region("simulation");
  lcaperf_.new_attribute("cycle",LCAP_INT);
  lcaperf_.new_attribute("level",LCAP_INT);

  parameters_ = new Parameters(parameter_file,monitor_);
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

Simulation::Simulation()
  : patch_loop_(0)
{ TRACE("Simulation()"); }

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

Simulation::Simulation (CkMigrateMessage *m)
  : patch_loop_(0)
{ TRACE("Simulation(CkMigrateMessage)"); }

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

  problem_->initialize_boundary(parameters_);
  problem_->initialize_initial (parameters_,group_process_);
  problem_->initialize_stopping(parameters_);
  problem_->initialize_timestep(parameters_);
  problem_->initialize_output  (parameters_,field_descr_,
				group_process_,factory());
  problem_->initialize_method  (parameters_);

  initialize_hierarchy_();

}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  DEBUG0;
  lcaperf_.stop("simulation");
  lcaperf_.print();
  lcaperf_.end();

  performance_simulation_->stop();
  performance_cycle_->stop();
}

//======================================================================

void Simulation::initialize_simulation_() throw()
{

  lcaperf_.attribute("cycle",&cycle_,LCAP_INT);
  lcaperf_.attribute("level",&level_,LCAP_INT);
  lcaperf_.start("simulation");
  performance_simulation_->start();
  performance_cycle_->start();

  //--------------------------------------------------
  // parameter: Mesh    : root_rank
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

  dimension_ = parameters_->value_integer("Mesh:root_rank",0);
  
  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be specified",
	  dimension_ != 0);
  
  ASSERT ("Simulation::initialize_simulation_()", 
	  "Parameter 'Mesh:root_rank' must be 1, 2, or 3",
	  (1 <= dimension_) && (dimension_ <= 3));

  cycle_ = parameters_->value_integer("Initial:cycle",0);
  time_  = parameters_->value_float  ("Initial:time",0);
  dt_ = 0;
}

//----------------------------------------------------------------------

void Simulation::initialize_monitor_() throw()
{
  //--------------------------------------------------
  // parameter: Monitor : debug
  //--------------------------------------------------

  bool debug = parameters_->value_logical("Monitor:debug",false);

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

  int i;
  for (i=0; i<parameters_->list_length("Field:fields"); i++) {
    field_descr_->insert_field
      (parameters_->list_value_string(i, "Field:fields"));
  }

  // Define default ghost zone depth for all fields, default value of 1

  //--------------------------------------------------
  // parameter: Field : ghosts
  //--------------------------------------------------

  int gx = 0;
  int gy = 0;
  int gz = 0;

  if (parameters_->type("Field:ghosts") == parameter_integer) {
    gx = gy = gz = parameters_->value_integer("Field:ghosts",0);
    if (dimension_ < 2) gy = 0;
    if (dimension_ < 3) gz = 0;
  } else if (parameters_->type("Field:ghosts") == parameter_list) {
    gx = parameters_->list_value_integer(0,"Field:ghosts",0);
    gy = parameters_->list_value_integer(1,"Field:ghosts",0);
    gz = parameters_->list_value_integer(2,"Field:ghosts",0);
  }

  for (i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_ghosts(i,gx,gy,gz);
  }

  // Set face dimensions to refresh

  //--------------------------------------------------
  // parameter: Field : refresh_faces
  // parameter: Field : refresh_edges
  // parameter: Field : refresh_corners
  //--------------------------------------------------

  // Refresh face ghost zones
  if (parameters_->type("Field:refresh_faces") == parameter_logical) {
    bool refresh_faces = 
      parameters_->value_logical ("Field:refresh:faces",true);
    field_descr_->set_refresh_face(2,refresh_faces);
  }

  // Refresh edge ghost zones
  if (parameters_->type("Field:refresh_edges") == parameter_logical) {
    bool refresh_edges = 
      parameters_->value_logical ("Field:refresh:edges",false);
    field_descr_->set_refresh_face(1,refresh_edges);
  }

  // Refresh corner ghost zones
  if (parameters_->type("Field:refresh_corners") == parameter_logical) {
    bool refresh_corners = 
      parameters_->value_logical ("Field:refresh:corners",false);
    field_descr_->set_refresh_face(0,refresh_corners);
  }
  
  //--------------------------------------------------
  // parameter: Field : precision
  //--------------------------------------------------

  std::string precision_str = 
    parameters_->value_string("Field:precision","default");

  precision_type precision = precision_unknown;

  if      (precision_str == "default")   precision = precision_default;
  else if (precision_str == "single")    precision = precision_single;
  else if (precision_str == "double")    precision = precision_double;
  else if (precision_str == "quadruple") precision = precision_quadruple;
  else {
    ERROR1 ("Simulation::initialize_data_descr_()", 
	    "Unknown precision %s",
	    precision_str.c_str());
  }

  for (i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_precision(i,precision);
  }

  //--------------------------------------------------
  // parameter: Field : alignment
  //--------------------------------------------------

  int alignment = parameters_->value_integer("Field:alignment",8);

  field_descr_->set_alignment (alignment);
  
  //--------------------------------------------------
  // parameter: Field : padding
  //--------------------------------------------------

  int padding = parameters_->value_integer("Field:padding",0);

  field_descr_->set_padding (padding);

  //--------------------------------------------------
  // parameter: Field : <field_name> : centering
  //--------------------------------------------------

  for (int i=0; i<field_descr_->field_count(); i++) {

    std::string field_name = field_descr_->field_name(i);

    std::string param_name = 
      std::string("Field:") + field_name + ":centering";
    
    if (parameters_->type(param_name) != parameter_unknown) {

      bool valid = parameters_->type(param_name) == parameter_list;
      valid = valid && parameters_->list_length(param_name) == dimension_;
      for (int i=0; i<dimension_; i++) {
	valid = valid && 
	  (parameters_->list_type(i,param_name) == parameter_logical);
      }
      
      ASSERT2 ("Simulation::initialize_data_descr_()",
	       "Parameter '%s' must be a list of logical values with length %d",
	       param_name.c_str(), dimension_, valid);

      int id_field = field_descr_->field_id(field_name);

      bool cx,cy,cz;

      cx = (dimension_ >= 1) ? 
	parameters_->list_value_logical(0,param_name,true) : true;
      cy = (dimension_ >= 2) ? 
	parameters_->list_value_logical(1,param_name,true) : true;
      cz = (dimension_ >= 3) ? 
	parameters_->list_value_logical(2,param_name,true) : true;

      field_descr_->set_centering (id_field, cx,cy,cz);

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
  hierarchy_ = factory()->create_hierarchy (dimension_,refinement);

  // Domain extents

  //--------------------------------------------------
  // parameter: Domain : lower
  // parameter: Domain : upper
  //--------------------------------------------------

  ASSERT ("Simulation::initialize_hierarchy_",
	  "Parameter Domain:lower list length must match Physics::dimension",
	  (parameters_->list_length("Domain:lower") == dimension_));

  ASSERT ("Simulation::initialize_hierarchy_",
	  "Parameter Domain:upper list length must match Physics::dimension",
	  (parameters_->list_length("Domain:upper") ==  dimension_));

  double lower[3];
  double upper[3];

  for (int i=0; i<3; i++) {
    lower[i] = parameters_->list_value_float(i, "Domain:lower", 0.0);
    upper[i] = parameters_->list_value_float(i, "Domain:upper", 0.0);
    ASSERT ("Simulation::initialize_hierarchy_",
	    "Domain:lower may not be greater than Domain:upper",
	    lower[i] <= upper[i]);
  }

  hierarchy_->set_lower(lower[0], lower[1], lower[2]);
  hierarchy_->set_upper(upper[0], upper[1], upper[2]);

  //----------------------------------------------------------------------
  // Create and initialize root Patch in Hierarchy
  //----------------------------------------------------------------------

  //--------------------------------------------------
  // parameter: Mesh : root_size
  // parameter: Mesh : root_blocks
  //--------------------------------------------------

  int root_size[3];

  root_size[0] = parameters_->list_value_integer(0,"Mesh:root_size",1);
  root_size[1] = parameters_->list_value_integer(1,"Mesh:root_size",1);
  root_size[2] = parameters_->list_value_integer(2,"Mesh:root_size",1);

  hierarchy_->set_root_size(root_size[0],root_size[1],root_size[2]);

  int root_blocks[3];

  root_blocks[0] = parameters_->list_value_integer(0,"Mesh:root_blocks",1);
  root_blocks[1] = parameters_->list_value_integer(1,"Mesh:root_blocks",1);
  root_blocks[2] = parameters_->list_value_integer(2,"Mesh:root_blocks",1);

#ifndef CONFIG_USE_CHARM
  ASSERT4 ("Simulation::initialize_hierarchy_",
	   "Product of Mesh:root_blocks = [%d %d %d] must equal MPI_Comm_size = %d",
	   root_blocks[0],root_blocks[1],root_blocks[2], group_process_->size(),
	   root_blocks[0]*root_blocks[1]*root_blocks[2]==group_process_->size());
#endif

  std::string type = parameters_->value_string("Initial:type","default");

  // Don't allocate blocks if reading data from files
  bool allocate_blocks = ! ( type == "file" || type == "restart" );

#ifdef CONFIG_USE_CHARM
  // Distributed patches in Charm: only allocate on root processor
  ++patch_loop_.stop();
  if (group_process()->is_root())
#endif
    {
      hierarchy_->create_root_patch
	(field_descr_,
	 root_size[0],root_size[1],root_size[2],
	 root_blocks[0],root_blocks[1],root_blocks[2],
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
  //  delete lcaperf_; lcaperf_ = 0;
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

void Simulation::update_cycle(int cycle, double time, double dt, double stop) 
{
  DEBUG4 ("Simulation::update_cycle cycle %d time %f dt %f stop %f",
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
  lcaperf_.print();
  lcaperf_.stop("simulation");


  lcaperf_.attribute("cycle",&cycle_,LCAP_INT);
  lcaperf_.start("simulation");
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
