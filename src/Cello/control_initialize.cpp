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

// #define DEBUG_INITIALIZE

#ifdef DEBUG_INITIALIZE
#  define TRACE_INITIALIZE CkPrintf ("%d %s:%d DEBUG_INITIALIZE\n",CkMyPe(),__FILE__,__LINE__);
#else
#  define TRACE_INITIALIZE /*  */
#endif

//----------------------------------------------------------------------

void Simulation::initialize() throw()
{
  TRACE_INITIALIZE;
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
  problem_->initialize_method(config_);
  problem_->initialize_initial(config_,parameters_);
  problem_->initialize_refine  (config_,parameters_);
  problem_->initialize_stopping(config_);
  problem_->initialize_output  (config_,factory());
  problem_->initialize_solver  (config_);
  problem_->initialize_prolong (config_);
  problem_->initialize_restrict (config_);

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
  TRACE_INITIALIZE;
  performance_->start_region(perf_initial);
  delete msg;
  
  initialize_block_array_();
}

//----------------------------------------------------------------------

void Simulation::r_initialize_hierarchy(CkReductionMsg * msg) 
{
  TRACE_INITIALIZE;
  performance_->start_region(perf_initial);
  delete msg;

  if (CkMyPe() == 0) {

    // --------------------------------------------------
    // ENTRY: #3 Simulation::r_initialize_hierarchy() -> Block::p_adapt_mesh()
    // ENTRY: Block Array if Simulation is_root()
    // --------------------------------------------------
    hierarchy()->block_array().p_initial_exit();
    // --------------------------------------------------
  }
  performance_->stop_region(perf_initial);
}

//======================================================================
// NEW INITIAL
//======================================================================

void Block::initial_enter_()
{
}

//----------------------------------------------------------------------

void  Block::initial_begin_()
{
}

//----------------------------------------------------------------------

void  Block::initial_next_()
{
}

//----------------------------------------------------------------------

void  Block::initial_continue_()
{
}

//----------------------------------------------------------------------

void  Block::initial_end_()
{
}

//----------------------------------------------------------------------

