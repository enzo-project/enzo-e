/** 
 *********************************************************************
 *
 * @file      test_hydro.cpp
 * @brief     Program implementing unit tests for hydrodynamics
 * @author    James Bordner
 * @date      Fri Mar  7 17:11:14 PST 2008
 *
 *********************************************************************
 */

#include "string.h"
#include "cello_hydro.h"
#include "test_ppm.h"
#include "performance_timer.hpp"

void initialize_hydro ();
void initialize_image ();
void initialize_implosion ();
void initialize_implosion2 ();


void print_usage(const char * name)
{
  printf ("Usage: %s <image|implosion|implosion3> [size] [cycles] [dump-frequency]\n",name);
}
int main(int argc, char ** argv)
{
  enum type_problem problem_type;
  int n;
  int cycles = 20000;
  int cycle_dump_frequency = 10;
  int block_size = 0;

  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }
  if (strcmp(argv[1],"image")==0) {
    problem_type = problem_image;
    printf ("image\n");
  } else if (strcmp(argv[1],"implosion")==0) {
    problem_type = problem_implosion;
  } else if (strcmp(argv[1],"implosion3")==0) {
    problem_type = problem_implosion3;
  } else {
    print_usage(argv[0]);
    exit(1);
  }
  if (argc>=3)  {
    n = atoi(argv[2]);
    if (n < 1 || 10000 < n) {
      int n_old = n;
      if (problem_type == problem_implosion3) {
	n = 32-6;
      } else {
	n = 400-6;
      }
      printf ("Illegal size %d: resetting to %d\n",n_old,n);
    }
  } else {
    if (problem_type == problem_implosion3) {
      n = 32-6;
    } else {
      n = 400-6;
    }
  }
  if (argc>=4)  {
    cycles = atoi(argv[3]);
    if (cycles < 1 || 10000000 < cycles) {
      printf ("Illegal cycles %d: resetting to 20000\n",n);
      cycles = 20000;
    }
  }
  if (argc>=5)  {
    cycle_dump_frequency = atoi(argv[4]);
    if (cycle_dump_frequency < 0) {
      printf ("Illegal cycle_dump_frequency %d: resetting to 10\n",n);
      cycle_dump_frequency = 10;
    }
  }
  if (argc>=6)  {
    block_size = atoi(argv[5]);
    if (block_size != 0 && (block_size < 4 || block_size > 256) ) {
      printf ("Illegal block_size %d: resetting to 0\n",block_size);
      block_size = 0;
    }
  }

  initialize_hydro ();

  switch (problem_type) {
  case problem_image:
    initialize_image(cycles);
    break;
  case problem_implosion:
    initialize_implosion(n,cycles);
    break;
  case problem_implosion3:
    initialize_implosion3(n,cycles);
    break;
  }

  float dt;

  int   cycle;
  float time;

  Timer timer;
  timer.start();
    
  for (cycle = 0, time = 0.0;
       (cycle < cycle_stop) && (time < time_stop);
       ++cycle, time += dt) {

    dt =  min(ComputeTimeStep(), time_stop - time);


    SetExternalBoundaryValues();

    if (cycle_dump_frequency && (cycle % cycle_dump_frequency) == 0) {
      printf ("cycle = %6d seconds = %5.0f sim-time = %6f dt = %6f\n",
	      cycle,timer.value(),time,dt);
      data_dump(problem_name[problem_type],cycle);
    }

    SolveHydroEquations(cycle, dt);

  }
  printf ("%d %d %g\n",n+6,cycles,timer.value());

  if (cycle_dump_frequency && (cycle % cycle_dump_frequency) == 0) {
    SetExternalBoundaryValues();
    data_dump(problem_name[problem_type],cycle);
  }
}

