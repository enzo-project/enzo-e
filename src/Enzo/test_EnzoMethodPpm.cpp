// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:04:03 PST 2008
/// @todo     Move initialization into Simulation / EnzoSimulation
/// @brief    Program implementing unit tests for the EnzoMethodPpm class

#include "test.hpp"

#include "enzo.hpp"
#include "test_EnzoMethodPpm.hpp"

#include PARALLEL_CHARM_INCLUDE(test_EnzoMethodPpm.decl.h)

PARALLEL_MAIN_BEGIN

{

  // Initialize parallelism

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  unit_init(parallel->rank(), parallel->size());

  // Initialize cross-cutting components

  Error     * error    = new Error;
  Monitor   * monitor  = new Monitor;

  // Create top-level Simulation object

  Simulation * simulation = new EnzoSimulation(error,monitor);

  FILE * fp = fopen ("input/test_EnzoMethodPpm.in","r");

  simulation->initialize(fp);

  fclose (fp);

  simulation->run();


  simulation->finalize();

  unit_finalize();

  delete error;
  delete monitor;
  delete simulation;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(test_EnzoMethodPpm.def.h)

