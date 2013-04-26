// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_SimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of SimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "simulation.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#include "simulation.hpp"

//----------------------------------------------------------------------

SimulationCharm::SimulationCharm
(
 const char         parameter_file[],
 int                n) throw ()
  : Simulation(parameter_file, n),
    block_sync_(0)
{
  TRACE("SimulationCharm::SimulationCharm");

  // derived class should call initialize()
}

//----------------------------------------------------------------------

SimulationCharm::~SimulationCharm() throw()
{
  TRACE("SimulationCharm::~SimulationCharm()");
}

//----------------------------------------------------------------------

void SimulationCharm::initialize() throw()
{
  TRACE("SimulationCharm::initialize()");
  Simulation::initialize();

  WARNING("SimulationCharm::initialize()",
	  "Calling StartLB for debugging load balancing()");

  s_initialize();
}

//----------------------------------------------------------------------

void SimulationCharm::s_initialize()
{
  TRACE("SimulationCharm::s_initialize()");
  if (group_process_->is_root()) run();
}

//----------------------------------------------------------------------

void SimulationCharm::run() throw()
{
  TRACE("SimulationCharm::run()");
  initial();
}

//----------------------------------------------------------------------

void SimulationCharm::p_refresh()
{
  TRACE("SimulationCharm::p_refresh");
  refresh();
};

void SimulationCharm::refresh()
{
  TRACE("SimulationCharm::refresh");
  if (hierarchy()->group_process()->is_root()) 
    hierarchy()->block_array()->p_refresh(); 
}

//----------------------------------------------------------------------

void SimulationCharm::c_compute()
{
  TRACE("SimulationCharm::c_compute()");
  if (stop_) {
    
    performance_.stop_region (id_cycle_);
    performance_write();

    proxy_main.p_exit(CkNumPes());

  } else {

    if (cycle_ > 0 ) performance_.stop_region (id_cycle_);

    performance_.start_region (id_cycle_);

    if (hierarchy()->group_process()->is_root()) 
      hierarchy()->block_array()->p_compute(cycle_,time_,dt_);
  }

}

//======================================================================

#endif /* CONFIG_USE_CHARM */
