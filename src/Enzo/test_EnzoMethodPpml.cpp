// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_EnzoMethodPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Apr  1 16:19:18 PDT 2010
/// @brief    Unit tests for the EnzoMethodPpml class

#include "test.hpp"

#include "enzo.hpp"

int main (int argc, char ** argv)
{

  // Initialize parallelism

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  unit_init(parallel->rank(), parallel->size());

  // Initialize cross-cutting components

  Monitor   * monitor  = new Monitor;

  // Create top-level Simulation object

  Simulation * simulation = new EnzoSimulation(monitor);

  FILE * fp = fopen ("input/test_EnzoMethodPpml.in","r");

  unit_assert(fp != 0);

  simulation->initialize(fp);

  fclose (fp);

  simulation->run();


  simulation->finalize();

  unit_finalize();

  delete monitor;
  delete simulation;

  PARALLEL_EXIT;
}
