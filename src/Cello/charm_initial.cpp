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

void SimulationCharm::initial ()
{
  problem()->initial_reset();
  problem()->initial_next(this);
}

//======================================================================

void Problem::initial_next(Simulation * simulation) throw()
{

  // find next initialization object

  Initial * initial;

  initial = this->initial(++index_initial_);

  Hierarchy *  hierarchy   = simulation->hierarchy();
  FieldDescr * field_descr = simulation->field_descr();

  if (initial != NULL) {

    if (initial->expects_blocks_allocated()) {

      DEBUG1 ("Start Initial(%d) A",index_initial_);

      TRACE ("DEBUG: Problem::initial_next calling CommBlock::p_initial()");
      if (hierarchy->group_process()->is_root()) {
	hierarchy->block_array()->p_initial();
      }

    } else {

      DEBUG1 ("Start Initial(%d) B",index_initial_);

      initial->enforce_block((CommBlock *)NULL, field_descr, hierarchy);

    }

  } else {

    ERROR("Problem::initial_next()",
	  "Multiple Initial objects not yet implemented");

    SimulationCharm * simulation_charm  = 
      dynamic_cast<SimulationCharm *> (proxy_simulation.ckLocalBranch());

    simulation_charm->refresh();

  }
}

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

  Initial * initial = simulation->problem()->initial();

  initial->enforce_block(this,field_descr, simulation->hierarchy());

  SimulationCharm * simulation_charm  = proxy_simulation.ckLocalBranch();
  simulation_charm->s_initial();

}

//----------------------------------------------------------------------

void SimulationCharm::s_initial()
{
  TRACE("ENTER SimulationCharm::s_initial()");
  TRACE2 ("block_sync: %d/%d",block_sync_.index(),block_sync_.stop());
  if (block_sync_.done()) {
    TRACE ("CONTINUE SimulationCharm::s_initial()");
    CkCallback callback (CkIndex_SimulationCharm::c_initial(), thisProxy);
    contribute(0,0,CkReduction::concat,callback);

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


