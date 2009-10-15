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

const char * file_root = "image";

int main(int argc, char **argv)
{

  initialize_hydro ();
  initialize_image();

  float dt;

  int   cycle;
  const int cycle_dump_frequency = 10;
  float time;

  for (cycle = 0, time = 0.0;
       (cycle < cycle_stop) && (time < time_stop);
       ++cycle, time += dt) {

    dt =  min (ComputeTimeStep(), time_stop - time);

    printf ("cycle = %6d time = %6f dt = %6f\n",cycle,time,dt);

    SetExternalBoundaryValues();

    if ((cycle % cycle_dump_frequency) == 0) {
      data_dump(file_root,cycle);
    }

    SolveHydroEquations(cycle, dt);

  }
  if ((cycle % cycle_dump_frequency) == 0) {
    SetExternalBoundaryValues();
    data_dump(file_root,cycle);
  }
}

