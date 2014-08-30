/***********************************************************************
/
/  GRID CLASS (SOLVE THE HYPERBOLIC EQUATIONS, SAVING SUBGRID FLUXES)
/
/  written by: Alexei Kritsuk
/  date:       May 2009
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Solve the hyperbilic equations with an array of solvers, saving the subgrid fluxes
//
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */

#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"


int grid::SolveHyperbolicEquations(int CycleNumber, int NumberOfSubgrids,
                              fluxes *SubgridFluxes[], int level)
{

  if (HydroMethod == PPML_Isothermal3D) {
    if (this->SolveMHDEquations(CycleNumber, NumberOfSubgrids, SubgridFluxes, level) == FAIL) {
      fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
      return FAIL;
    }
  }
  else
    if (this->SolveHydroEquations(CycleNumber, NumberOfSubgrids, SubgridFluxes, level) == FAIL) {
      fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
      return FAIL;
    }

  return SUCCESS; 
  
}
