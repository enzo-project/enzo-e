// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      test_ppm.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @todo      Replace image_dump() with Monitor::image()
/// @date      Fri Mar  7 17:11:14 PST 2008
/// @brief     Program implementing unit tests for hydrodynamics

#include "test.hpp"

#include "enzo.hpp"

const int num_problems = 4;
const int num_ghosts   = 3;

enum problem_ppm_enum {  
  problem_ppm_unknown, 
  problem_ppm_image,  
  problem_ppm_implosion,  
  problem_ppm_implosion3
};

const char * problem_name[] = {
  "",
  "ppm-image",
  "ppm-implosion",
  "ppm-implosion3"
};

const int problem_size [] = {
  0,
  512,
  400,
  32 
};

const int problem_cycles [] = {
  0,
  10000,
  10000,
  10000 
};

#include PARALLEL_CHARM_INCLUDE(test_ppm.decl.h)

//----------------------------------------------------------------------

void print_usage(const char * name)
{
  PARALLEL_PRINTF ("Usage: %s <ppm-image|ppm-implosion|ppm-implosion3> [size] [cycles] [dump-frequency]\n",name);
  exit(1);
}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN

{

  PARALLEL_INIT;

  unit_init();

  // Create global objects

  Error   * error   = new Error;
  Monitor * monitor = new Monitor;

  // Initialize monitor

  // Check command line arguments

  if (PARALLEL_ARGC < 2) {
    print_usage(PARALLEL_ARGV[0]);
  }

  int argi = 0;

  // Parse command line arguments

  enum problem_ppm_enum problem = problem_ppm_unknown;

  if (PARALLEL_ARGC > ++argi) {
    for (int i=0; i<num_problems; i++) {
      if (strcmp(PARALLEL_ARGV[argi],problem_name[i]) == 0) {
	problem = problem_ppm_enum(i);
      }
    }
  }
  if (problem == 0) print_usage(PARALLEL_ARGV[0]);

  int size = problem_size[problem];
  if (PARALLEL_ARGC > ++argi) size = atoi(PARALLEL_ARGV[argi]);
  int cycle_stop = problem_cycles[problem];
  if (PARALLEL_ARGC > ++argi) cycle_stop = atoi(PARALLEL_ARGV[argi]);
  int dump_frequency = 10;
  if (PARALLEL_ARGC > ++argi) dump_frequency = atoi(PARALLEL_ARGV[argi]);

  PARALLEL_PRINTF ("problem = %s  size = %d  cycles = %d  dump_frequency = %d\n",
	  problem_name[problem], size, cycle_stop, dump_frequency);

  // Initialize for generic hydrodynamics

  EnzoDescr enzo;
  enzo.initialize_hydro ();

  // Initialize for specific problem type

  float  time_stop  = 2.5;

  switch (problem) {
  case problem_ppm_image:
    enzo.initialize_image();
    break;
  case problem_ppm_implosion:
    enzo.initialize_implosion(size);
    break;
  case problem_ppm_implosion3:
    enzo.initialize_implosion3(size);
    break;
  default:
    print_usage(PARALLEL_ARGV[0]);
  }

  float dt;

  int   cycle;
  float time;

  Papi papi;
  Timer timer;
  timer.start();
  papi.start();

  double lower = 0.125*size;
  double upper =   1.0*size;

  unit_class ("enzo");
  unit_func ("SolveHydroEquations");

  for (cycle = 0, time = 0.0;
       (cycle < cycle_stop) && (time < time_stop);
       ++cycle, time += dt) {

    enzo.SetExternalBoundaryValues();
    dt =  MIN(enzo.ComputeTimeStep(), time_stop - time);

    if (dump_frequency && (cycle % dump_frequency) == 0) {
      PARALLEL_PRINTF ("cycle = %6d seconds = %5.0f sim-time = %10g dt = %10g\n",
	      cycle,timer.value(),time,dt);
      fflush(stdout);
      enzo.image_dump(problem_name[problem],cycle,lower,upper,monitor);
    }

    enzo.SolveHydroEquations(NULL,cycle, dt);

  }
  unit_assert(cycle >=cycle_stop || time >= time_stop);

  timer.stop();
  papi.stop();

  PARALLEL_PRINTF ("(size+6,cycles,time) = %d %d %g\n",
		   size+6,cycle_stop,timer.value());
  fflush(stdout);

  if (dump_frequency && (cycle % dump_frequency) == 0) {
    enzo.SetExternalBoundaryValues();
    enzo.image_dump(problem_name[problem],cycle,lower,upper,monitor);
  }

  printf ("Time real = %f\n",papi.time_real());
  printf ("Time proc = %f\n",papi.time_proc());
  printf ("Flop count = %lld\n",papi.flop_count());
  printf ("GFlop rate = %f\n",papi.flop_rate()*1e-9);

  delete monitor;

  unit_finalize();

  PARALLEL_EXIT;

}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_ppm.def.h)
