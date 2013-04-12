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
    block_loop_(0)
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

#ifdef REMOVE_PATCH

  s_initialize();

#endif

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
#ifdef REMOVE_PATCH
  hierarchy()->block_array()->p_refresh(); 
#else /* REMOVE_PATCH */
  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while (( patch = ++it_patch )) {
    CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
    proxy_patch->p_refresh();
  }
#endif /* REMOVE_PATCH */

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

#ifdef REMOVE_PATCH
    hierarchy()->block_array()->p_compute(cycle_,time_,dt_);
#else /* REMOVE_PATCH */
    ItPatch it_patch(hierarchy_);
    Patch * patch;
    while (( patch = ++it_patch )) {
      CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
      proxy_patch->p_compute(cycle_, time_, dt_);
    }
#endif /* REMOVE_PATCH */

  }

}

//======================================================================

#endif /* CONFIG_USE_CHARM */
