// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Simulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @todo     Add create_simulation
/// @brief    Program implementing unit tests for the Simulation class
 
#include "test.hpp"

#include "simulation.hpp"
#include "enzo.hpp" /* Required for EnzoSimulationSerial */

#include PARALLEL_CHARM_INCLUDE(test_Simulation.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  GroupProcess * group_process = GroupProcess::create();

  unit_init();

  // Create the simulation

  unit_func("Simulation","Simulation");

  // NOTE: Need concrete EnzoSimulationSerial class since Simulation is 
  //       an abstract base class

  FILE * fp = fopen ("input/test_Simulation.in","r");
  Parameters * parameters = new Parameters;
  parameters -> read(fp);
  fclose (fp);

  Simulation * simulation = new EnzoSimulationMpi(parameters,group_process);

  unit_assert(simulation != 0);

  // Initialize the simulation

  unit_func("Simulation","initialize");

  simulation->initialize();

  unit_assert(true);

  // Test accessor functions

  unit_func("Simulation","mesh");
  unit_assert (simulation->mesh() != NULL);

  unit_func("Simulation","stopping");
  unit_assert (simulation->stopping() != NULL);

  unit_func("Simulation","timestep");
  unit_assert (simulation->timestep() != NULL);

  unit_func("Simulation","initial");
  unit_assert (simulation->initial() != NULL);

  unit_func("Simulation","boundary");
  unit_assert (simulation->boundary() != NULL);

  unit_func("Simulation","num_method");
  unit_assert (simulation->num_method() > 0);

  unit_func("Simulation","method");
  for (int i=0; i<simulation->num_method(); i++) {
    unit_assert (simulation->method(i) != NULL);
  }

  // Run the simulation

  unit_func("Simulation","run");
  simulation->run();
  unit_assert(unit_incomplete);
  
  // Load the simulation

  unit_func("Simulation","read");
  simulation->read();

  // Save the simulation

  unit_func("Simulation","write");
  simulation->write();
  unit_assert(unit_incomplete);

  // Finalize the simulation

  unit_func("Simulation","finalize");
  simulation->finalize();
  unit_assert(unit_incomplete);

  delete simulation;

  unit_finalize();


  PARALLEL_EXIT;

}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Simulation.def.h)
