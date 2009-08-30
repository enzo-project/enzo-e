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

main()
{

  initialize_hydro ();
  initialize_implosion();

  float dt;

  int   cycle;
  float time;

  for (cycle = 0, time = 0.0;
       (cycle < cycle_stop) && (time < time_stop);
       ++cycle, time += dt) {

    dt =  min (ComputeTimeStep(), time_stop - time);

    printf ("%s:%d\n",__FILE__,__LINE__);
    printf ("   dt = %g\n",dt);
    printf (" sum (density) = %g\n",sum_grid(field_density));
    printf (" sum (total_energy) = %g\n",sum_grid(field_total_energy));
    printf (" sum (velocity_x) = %g\n",sum_grid(field_velocity_x));
    printf (" sum (velocity_y) = %g\n",sum_grid(field_velocity_y));


    SolveHydroEquations(cycle, dt);

  }
}
