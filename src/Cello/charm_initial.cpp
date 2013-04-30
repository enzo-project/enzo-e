// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_charm_initial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-04-03
/// @brief    Functions implementing initialization functions requiring
///             CHARM++
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
  // Apply the initial conditions 

  apply_initial_();

  CkStartQD (CkCallback(CkIndex_CommBlock::p_phase_adapt(),thisProxy[thisIndex]));

}

//----------------------------------------------------------------------

void CommBlock::apply_initial_() throw ()
{

  SimulationCharm * simulation  = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  // Set the CommBlock cycle and time to match Simulation's

  TRACE("CommBlock::p_initial Setting time");
  set_cycle(simulation->cycle());
  set_time (simulation->time());
  set_dt   (simulation->dt());

  // Allocate block data

  block()->allocate(field_descr);

  // Perform any additional initialization for derived class 

  initialize ();

  // Apply initial conditions

  index_initial_ = 0;
  while (Initial * initial = simulation->problem()->initial(index_initial_++)) {
    initial->enforce_block(this,field_descr, simulation->hierarchy());
  }
}

//----------------------------------------------------------------------

void CommBlock::p_read (int index_initial)
{
  INCOMPLETE("CommBlock::p_read");
}

//======================================================================

#endif /* CONFIG_USE_CHARM */


