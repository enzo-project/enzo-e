// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoStopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-08-28
/// @brief    Implementation of EnzoStopping class for Enzo-specific stopping criteria

#include "charm_simulation.hpp"
#include "enzo.hpp"

//======================================================================

bool EnzoStopping::complete (int    curr_cycle,
			   double curr_time) const throw()
{
  bool stop = Stopping::complete(curr_cycle,curr_time);

  // check for stopping on redshift

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  EnzoUnits * units = (EnzoUnits * )simulation->problem()->units();
  EnzoPhysicsCosmology * cosmology = units->cosmology();

  if (cosmology) {
    cosmology->set_current_time(curr_time);
    stop = (stop || (cosmology->current_redshift() <= stop_redshift_));
  }
    
  return stop;
}
