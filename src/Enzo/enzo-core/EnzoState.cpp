// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoState
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2024-06-14
/// @brief    State object including Enzo state, e.g. redshift

#include "Enzo/cosmology/cosmology.hpp"
#include "Enzo/enzo.hpp"

void EnzoState::set_time (double time)
{
  time_ = time;
  Simulation * simulation = cello::simulation();
  EnzoUnits * units = (EnzoUnits * )simulation->problem()->units();
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  if (cosmology) {
    cosmology->set_current_time(time);
    redshift_ = cosmology->current_redshift();
  }
}
