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
  : Simulation(parameter_file, n)
{
  // derived class should call initialize()
}

//----------------------------------------------------------------------

SimulationCharm::~SimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void SimulationCharm::initialize() throw()
{
  Simulation::initialize();

  WARNING("SimulationCharm::initialize()",
	  "Calling StartLB for debugging load balancing()");

}

//----------------------------------------------------------------------

void SimulationCharm::s_initialize()
{
  if (patch_loop_.done()) run();
}

//----------------------------------------------------------------------

void SimulationCharm::run() throw()
{
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
  ItPatch it_patch(hierarchy_);
  Patch * patch;

  while (( patch = ++it_patch )) {
    CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
    proxy_patch->p_refresh();
  }
}

//----------------------------------------------------------------------

void SimulationCharm::c_compute()
{
  if (stop_) {
    
    performance_.stop_region (id_cycle_);

    proxy_main.p_exit(CkNumPes());

  } else {

    if (cycle_ > 0 ) performance_.stop_region (id_cycle_);

    performance_.start_region (id_cycle_);

    ItPatch it_patch(hierarchy_);
    Patch * patch;
    while (( patch = ++it_patch )) {
      CProxy_Patch * proxy_patch = (CProxy_Patch *)patch;
      proxy_patch->p_compute(cycle_, time_, dt_);
    }
  }

}

//======================================================================

#endif /* CONFIG_USE_CHARM */
