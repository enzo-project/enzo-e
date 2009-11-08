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
 
#include "cello_hydro.h"

#include "performance.hpp"

const char * file_root = "ppml";

int main(int argc, char * argv[])
{

  int n      = 64;
  int cycles = 20000;
  int cycle_dump_frequency = 10;
  if (argc>=2)  {
    n = atoi(argv[1]);
    if (n < 1 || 200 < n) {
      printf ("Illegal size %d: resetting to 64\n",n);
      n = 64;
    }
  }
  if (argc>=3)  {
    cycles = atoi(argv[2]);
    if (cycles < 1 || 10000000 < cycles) {
      printf ("Illegal cycles %d: resetting to 20000\n",n);
      cycles = 20000;
    }
  }
  if (argc>=4)  {
    cycle_dump_frequency = atoi(argv[3]);
    if (cycle_dump_frequency < 0) {
      printf ("Illegal cycle_dump_frequency %d: resetting to 10\n",n);
      cycle_dump_frequency = 10;
    }
  }

  printf ("%d %d %d\n",n,cycles,cycle_dump_frequency);

  initialize_hydro ();
  initialize_ppml(n,cycles);

  float dt;

  int   cycle;
  float time;

  Timer timer;
  timer.start();
  for (cycle = 0, time = 0.0;
       (cycle < cycle_stop) && (time < time_stop);
       ++cycle, time += dt) {

    dt =  min (ComputeTimeStep(), time_stop - time);

    printf ("cycle = %6d time = %6f dt = %6f\n",cycle,time,dt);

    SetExternalBoundaryValues();

    if ((cycle % cycle_dump_frequency) == 0) {
      data_dump(file_root,cycle);
    }

    SolveMHDEquations(cycle, dt);

  }
  if ((cycle % cycle_dump_frequency) == 0) {
    SetExternalBoundaryValues();
    data_dump(file_root,cycle);
  }
}

