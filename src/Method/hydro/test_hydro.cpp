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



main()
{
  unit_class ("MethodPpm");
  unit_open();
  grid g;

  initialize_hydro ();

  int      cycle = 0;
  int      num_subgrids = 0;
  fluxes * subgrid_fluxes[1] = {NULL};
  int      level = 0;

  g.SolveHydroEquations(cycle, num_subgrids, subgrid_fluxes, level);
}
