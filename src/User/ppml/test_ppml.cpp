// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file      test_ppml.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @todo      Replace image_dump() with Monitor::image()
/// @date      Fri Mar  7 17:11:14 PST 2008
/// @brief     Program implementing unit tests for hydrodynamics

#include <mpi.h>

#include "cello.hpp" 
#include "cello_hydro.h"

#include "enzo.hpp"
#include "string.h"
#include "test_ppml.h"
#include "parallel.hpp"
#include "performance.hpp"
#include "monitor.hpp"

//----------------------------------------------------------------------

void print_usage(const char * name)
{
  printf ("usage: %s <ppml-blast|ppml-implosion3> [size] [cycles] [dump-frequency]\n",
	  name);
  exit(1);
}

//----------------------------------------------------------------------

int main(int argc, char * argv[])
{

  Monitor * monitor = new Monitor;

  // Check command line arguments

  if (argc < 2) {
    print_usage(argv[0]);
  }

  int argi = 0;

  // Parse command line arguments

  enum problem_ppml_type problem = problem_ppml_unknown;

  if (argc > ++argi) {
    for (int i=0; i<num_problems; i++) {
      if (strcmp(argv[argi],problem_name[i]) == 0) {
	problem = problem_ppml_type(i);
      }
    }
  }
  if (problem == 0) print_usage(argv[0]);

  int size = problem_size[problem];
  if (argc > ++argi) size = atoi(argv[argi]);
  int cycle_stop = problem_cycles[problem];
  if (argc > ++argi) cycle_stop = atoi(argv[argi]);
  int dump_frequency = 10;
  if (argc > ++argi) dump_frequency = atoi(argv[argi]);

  printf ("problem = %s  size = %d  cycles = %d  dump_frequency = %d\n",
	  problem_name[problem], size, cycle_stop, dump_frequency);

  // Initialize for generic hydrodynamics

  EnzoDescr enzo;

  enzo.initialize_hydro ();

  // Initialize for specific problem type

  float  time_stop  = 2.5;

  switch (problem) {
  case problem_ppml_blast:
    enzo.initialize_ppml(size);
    break;
  case problem_ppml_implosion3:
    enzo.initialize_ppml_implosion3(size);
    break;
  default:
    print_usage(argv[0]);
  }

  float dt;

  int   cycle;
  float time;

  Timer timer;
  timer.start();

  double lower = 0.125*size;
  double upper =   1.0*size;

  for (cycle = 0, time = 0.0;
       (cycle < cycle_stop) && (time < time_stop);
       ++cycle, time += dt) {

    dt =  MIN(enzo.ComputeTimeStep(), time_stop - time);

    enzo.SetExternalBoundaryValues();

    if (dump_frequency && (cycle % dump_frequency) == 0) {
      printf ("cycle = %6d seconds = %5.0f sim-time = %6f dt = %6f\n",
	      cycle,timer.value(),time,dt);
      fflush(stdout);
      enzo.image_dump(problem_name[problem],cycle,lower,upper,monitor);
    }

    enzo.SolveMHDEquations(cycle, dt);

  }

  printf ("%d %d %g\n",size+6,cycle_stop,timer.value());
  fflush(stdout);

  if (dump_frequency && (cycle % dump_frequency) == 0) {
    enzo.SetExternalBoundaryValues();
    enzo.image_dump(problem_name[problem],cycle,lower,upper,monitor);
  }
}

