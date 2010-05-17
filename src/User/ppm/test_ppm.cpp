// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file      test_hydro.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @todo      Replace image_dump() with Monitor::image()
/// @date      Fri Mar  7 17:11:14 PST 2008
/// @brief     Program implementing unit tests for hydrodynamics

#include "cello.h"
#include "cello_hydro.h"
#include "test_ppm.h"

#include "string.h"
#include "parallel.hpp"
#include "performance.hpp"

//----------------------------------------------------------------------

void print_usage(const char * name)
{
  printf ("Usage: %s <image|implosion|implosion3> [size] [cycles] [dump-frequency]\n",name);
  exit(1);
}

//----------------------------------------------------------------------

int main(int argc, char ** argv)
{

  // Initialize parallelism

  Parallel::instance()->initialize(&argc,&argv);

  // Check command line arguments

  enum type_problem problem_type = problem_unknown;

  if (argc < 2) {
    print_usage(argv[0]);
  }

  int argi = 0;

  // Parse command line arguments

  if (argc > ++argi) {
    for (int i=0; i<num_problems; i++) {
      if (strcmp(argv[argi],problem_name[i]) == 0) {
	problem_type = type_problem(i);
      }
    }
  }
  if (problem_type == 0) print_usage(argv[0]);

  int size = problem_size[problem_type];
  if (argc > ++argi) size = atoi(argv[argi]);
  int cycles = problem_cycles[problem_type];
  if (argc > ++argi) cycles = atoi(argv[argi]);
  int dump_frequency = 10;
  if (argc > ++argi) dump_frequency = atoi(argv[argi]);

  printf ("problem = %s  size = %d  cycles = %d  dump_frequency = %d\n",
	  problem_name[problem_type], size, cycles, dump_frequency);

  // Initialize for generic hydrodynamics

  initialize_hydro ();

  // Initialize for specific problem type

  switch (problem_type) {
  case problem_image:
    initialize_image(cycles);
    break;
  case problem_implosion:
    initialize_implosion(size,cycles);
    break;
  case problem_implosion3:
    initialize_implosion3(size,cycles);
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

    dt =  min(ComputeTimeStep(), time_stop - time);

    SetExternalBoundaryValues();

    if (dump_frequency && (cycle % dump_frequency) == 0) {
      printf ("cycle = %6d seconds = %5.0f sim-time = %6f dt = %6f\n",
	      cycle,timer.value(),time,dt);
      fflush(stdout);
      image_dump(problem_name[problem_type],cycle,lower,upper);
    }

    SolveHydroEquations(cycle, dt);

  }

  printf ("%d %d %g\n",size+6,cycles,timer.value());
  fflush(stdout);

  if (dump_frequency && (cycle % dump_frequency) == 0) {
    SetExternalBoundaryValues();
    image_dump(problem_name[problem_type],cycle,lower,upper);
  }
}

