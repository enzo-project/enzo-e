// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_initial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-04-03
/// @brief    Functions implementing initialization functions requiring
/// CHARM++
///
/// This file contains member functions for various CHARM++ chares and
/// classes used for calling Initial objects in a CHARM++ simulation.
/// Functions are listed in roughly the order of flow-of-control.
///
///    INITIAL
///
///    SimulationCharm::initial()
///       problem.initial_reset()
///       problem.initial_next()
///
///    Problem::initial_next()
///       if (initial)
///          if (blocks_allocated)
///             if (is_root())
///                block_array.p_initial()
///          else
///             inital.enforce_block(0)
///       else [ NOT CALLED ]
///           >>>>> simulation_charm.refresh() >>>>>
///     
///    CommBlock::p_initial()
///       block()->allocate()
///       set (cycle,time,dt)
///       initialize()
///       initial->enforce_block(this)
///       simulation_charm.s_initial()
///
///    SimulationCharm::s_initial()
///       if (block_sync_.done())
///          contribute(SimulationCharm::c_initial())
///
///    SimulationCharm::c_initial()
///       delete parameters_
///       if (is_root()) 
///          >>>>> block_array->p_adapt(0) >>>>> 
///             
///    ----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"


//----------------------------------------------------------------------

void CommBlock::p_initial()
{
  TRACE("CommBlock::p_initial()");
  Simulation * simulation  = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  // Initialize the block

  block()->allocate(field_descr);

  // Set the CommBlock cycle and time to match Simulation's

  TRACE("CommBlock::p_initial Setting time");
  set_cycle(simulation->cycle());
  set_time (simulation->time());
  set_dt   (simulation->dt());

  // Perform any additional initialization for derived class 

  initialize ();

  // Apply the initial conditions 

  apply_initial_();

}

//----------------------------------------------------------------------

void CommBlock::apply_initial_() throw ()
{
  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation_charm->field_descr();
  index_initial_ = 0;
  while (Initial * initial = simulation_charm->problem()->initial(index_initial_++)) {
    initial->enforce_block(this,field_descr, simulation_charm->hierarchy());
  }
  proxy_simulation.s_initial();

  CkStartQD (CkCallback(CkIndex_CommBlock::p_phase_adapt(),thisProxy[thisIndex]));

}

//----------------------------------------------------------------------

void SimulationCharm::s_initial()
{
  TRACE("ENTER SimulationCharm::s_initial()");
  TRACE2 ("block_sync: %d/%d",block_sync_.index(),block_sync_.stop());
  if (block_sync_.done()) {
    TRACE ("CONTINUE SimulationCharm::s_initial()");
    // DOES NOT WORK HERE: race condition can lead to 
    //    parameters_ being accessed after being deleted
    //    delete parameters_;
    //    parameters_ = 0;

  }
}
//----------------------------------------------------------------------

void SimulationCharm::c_initial()
{
  TRACE("SimulationCharm::c_initial()");
  delete parameters_;
  parameters_ = 0;
  if (hierarchy()->group_process()->is_root()) 
    hierarchy()->block_array()->p_adapt(0); 
}

//----------------------------------------------------------------------

void CommBlock::p_read (int index_initial)
{
  INCOMPLETE("CommBlock::p_read");
}

//======================================================================

#endif /* CONFIG_USE_CHARM */


