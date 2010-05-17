// $Id: test_simulation.cpp 1302 2010-03-17 00:16:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_simulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Simulation class
 
#include <stdio.h>
#include <string>

#include "cello.h"
#include "test.hpp"
#include "simulation.hpp"

int main(int argc, char ** argv)
{
  unit_init();
  unit_class ("Simulation");
  unit_func("Simulation");
  Simulation simulation;
  unit_assert(false);
  unit_finalize();
}
