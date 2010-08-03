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
#include "simulation.hpp"

int main(int argc, char ** argv)
{
  Mpi::init(&argc,&argv);

  // Warning: unused
  GroupProcess * parallel = new GroupProcessMpi;

  unit_init();
  unit_class ("Simulation");
  unit_func("Simulation");
  Global * global = new Global;

  // Read parameter file

  FILE * fp = fopen("input/implosion.in","r");
  global->parameters()->read(fp);

  // Create simulation object

  Simulation simulation(global);

  // Initialize data fields

  simulation.initialize();

  double stop_time  = 1.0;
  int    stop_cycle = 1000;

  // Advance the simulation to the given time or cycle limit

  simulation.advance(stop_time,stop_cycle);
  
  unit_assert(false);
  unit_finalize();

  Mpi::finalize();
}
