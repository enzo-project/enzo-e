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

  // Read parameter file

  FILE * fp = fopen("input/implosion.in","r");
  global->parameters()->read(fp);

  // Create simulation object

  Global * global = new Global;

  // Create the simulation

  unit_func("Simulation");
  Simulation simulation(global);
  unit_assert(true);

  // Initialize the simulation

  unit_func("initialize");
  simulation.initialize();
  unit_assert(true);

  // Run the simulation

  unit_func("execute");
  simulation.execute();
  unit_assert(false);
  
  // Load the simulation

  unit_func("load");
  simulation.load();

  // Save the simulation

  unit_func("save");
  simulation.save();
  unit_assert(false);

  unit_finalize();


  PARALLEL_EXIT;

}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Simulation.def.h)
