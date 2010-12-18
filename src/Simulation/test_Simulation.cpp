// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Simulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Simulation class
 
#include <stdio.h>
#include <string>

#include "cello.hpp"

#include "test.hpp"
#include "enzo.hpp"
#include "simulation.hpp"
#include "parallel.def"

#include PARALLEL_CHARM_INCLUDE(test_Simulation.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class ("Simulation");

  // Create simulation object

  Global * global = new Global;

  // Create the simulation

  unit_func("Simulation");

  TRACE_MESSAGE;
  // NOTE: Need concrete EnzoSimulation class since Simulation is 
  //       an abstract base class

  Simulation * simulation = new EnzoSimulation (global);

  unit_assert(true);

  // Initialize the simulation

  unit_func("initialize");
  TRACE_MESSAGE;
  simulation->initialize("input/implosion.in");
  TRACE_MESSAGE;
  unit_assert(true);

  // Run the simulation

  unit_func("run");
  simulation->run();
  unit_assert(false);
  
  // Load the simulation

  unit_func("read");
  simulation->read();

  // Save the simulation

  unit_func("write");
  simulation->write();
  unit_assert(false);

  // Finalize the simulation

  unit_func("finalize");
  simulation->finalize();
  unit_assert(false);

  delete simulation;

  unit_finalize();


  PARALLEL_EXIT;

}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Simulation.def.h)
