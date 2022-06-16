// See LICENSE_CELLO file for license and copyright information

/// @file     control_initialize.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_INITIAL

#ifdef DEBUG_INITIAL
#   define TRACE_INITIAL(MSG,BLOCK)                     \
  CkPrintf ("TRACE_CONTROL_INITIAL %s %s\n",                    \
            BLOCK->name().c_str(),MSG); fflush(stdout);

#   define TRACE_INITIAL_SIM(MSG)                       \
  CkPrintf ("TRACE_CONTROL_INITIAL %s\n",MSG); fflush(stdout);
#else
#   define TRACE_INITIAL(MSG,BLOCK) /* ... */
#   define TRACE_INITIAL_SIM(MSG) /* ... */
#endif

//----------------------------------------------------------------------

void Simulation::initialize() throw()
{
  TRACE_INITIAL_SIM("Simulation::initialize()");
  set_phase(phase_initial);

  initialize_config_();

  initialize_monitor_();
  initialize_memory_();
  initialize_performance_();
  initialize_simulation_();
  initialize_balance_();

  initialize_data_descr_();

  problem_->initialize_units (config_);
  problem_->initialize_physics (config_,parameters_);
  problem_->initialize_boundary(config_,parameters_);
  problem_->initialize_prolong (config_);
  problem_->initialize_restrict (config_);
  problem_->initialize_initial(config_,parameters_);
  problem_->initialize_method  (config_,factory());
  problem_->initialize_solver  (config_);
  problem_->initialize_refine  (config_,parameters_);
  problem_->initialize_stopping(config_);
  problem_->initialize_output  (config_,factory());

  cello::finalize_fields();

  initialize_hierarchy_();

  // initialize_block_array() is called in charm_initialize
  // using QD to ensure that initialize_hierarchy() is called
  // on all processors before Blocks are created

  // Create the Block chare array
  CProxy_Block block_array;
  if (CkMyPe() == 0) {
    bool allocate_data = true;
    block_array = hierarchy_->new_block_proxy (allocate_data);
    thisProxy.p_set_block_array(block_array);
  }

  CkCallback callback 
    (CkIndex_Simulation::r_initialize_block_array(NULL), thisProxy);

  // --------------------------------------------------
#ifdef TRACE_CONTRIBUTE  
  CkPrintf ("%s:%d DEBUG_CONTRIBUTE r_initialize_block_array()\n",__FILE__,__LINE__); fflush(stdout);
#endif  
  contribute(0,0,CkReduction::concat,callback);
  // --------------------------------------------------
}

//----------------------------------------------------------------------

void Simulation::r_initialize_block_array(CkReductionMsg * msg) 
{
  TRACE_INITIAL_SIM("Simulation::r_initialize_block_array_()");
  performance_->start_region(perf_initial);
  delete msg;
  
  initialize_block_array_();
}

//======================================================================
// NEW INITIAL
//======================================================================

void  Block::initial_new_begin_(int level)
{
  TRACE_INITIAL("initial_new_begin_",this);
  index_initial_ = 0;

  CkCallback callback 
    (CkIndex_Block::r_initial_new_next(nullptr), thisProxy);

  contribute(0,0,CkReduction::concat,callback);
}

//----------------------------------------------------------------------

void  Block::initial_new_next_()
{
  TRACE_INITIAL("initial_new_next_()",this);
  Initial * initial = cello::problem()->initial(index_initial_);
  if (initial) {
    initial->enforce_block(this,nullptr);
  } else {
    bool is_first_cycle = (cycle_ == cello::config()->initial_cycle);
    if (is_first_cycle && level() <= 0) {
      initial_exit_();
    }
  }
}

//----------------------------------------------------------------------

void  Block::initial_done()
{
  if (cello::config()->initial_new) {

    TRACE_INITIAL("initial_new_done",this);

    // barrier before continuing
    CkCallback callback (CkIndex_Block::r_initial_new_continue(nullptr), thisProxy);

    contribute(0,0,CkReduction::concat,callback);
  }
}

//----------------------------------------------------------------------

void  Block::initial_new_continue_()
{
  TRACE_INITIAL("initial_new_continue",this);
  index_initial_++;
  initial_new_next_();
}

//----------------------------------------------------------------------


