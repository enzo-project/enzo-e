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
 
#include <stdio.h>

#include "test.hpp"

#include "cello_hydro.h"

bool is_converged 
( 
 double time,
 double time_stop,
 int    cycle,
 int    cycle_stop
) 
{
  return (time >= time_stop || cycle >= cycle_stop);
}

main()
{

  initialize_hydro ();
  initialize_implosion();

  float dt;

  int   cycle = 0;
  float time = 0.0;

  grid g;

  for (cycle = 0, time = 0.0;
       ! is_converged (time,time_stop,cycle,cycle_stop);
       ++cycle, time += dt) {

    dt = g.ComputeTimeStep();
    dt = min (dt, time_stop - time);

    g.SolveHydroEquations(cycle, dt);

  }
}
