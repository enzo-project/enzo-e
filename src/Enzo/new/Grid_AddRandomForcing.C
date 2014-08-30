/***********************************************************************
/
/  GRID CLASS (ADD RANDOM FORCING FIELDS TO VELOCITIES)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::AddRandomForcing(float * norm, float * bulkMomentum, float dtTopGrid)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "GARF: Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
 
  /* error check. */
 
  if (RandomForcingField[0] == NULL)
    ERROR_MESSAGE;
 
  if (RandomForcingField[0][0] == 0.0)
    ERROR_MESSAGE;
 
  if (dtTopGrid == 0.0)
    ERROR_MESSAGE;
 
 
  /* compute the field size */
 
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* update total energy first. */
 
 
  float levelNorm = (*norm)*dtFixed/dtTopGrid;
  if (debug & levelNorm < 0.0)
    WARNING_MESSAGE;
 
  /* do not do the update if using ZEUS */
 
  if (HydroMethod != Zeus_Hydro)
    for (int i = 0; i < size; i++)
      for (int dim = 0; dim < GridRank; dim++)
	BaryonField[TENum][i] +=
	  (BaryonField[Vel1Num+dim][i]*
	   //	  RandomForcingField[dim][i])*levelNorm;
  	  (RandomForcingField[dim][i]-bulkMomentum[dim]))*levelNorm;
 
  /* add velocity perturbation to the velocity fields;
     keep the center-of-mass velocity zero. */
 
  if (debug && GridRank == 3)
    fprintf(stderr, "GARF: Bulk Momentum input %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", 
	    bulkMomentum[0], bulkMomentum[1], bulkMomentum[2], levelNorm);
      

  for (int dim = 0; dim < GridRank; dim++)
    for (int i = 0; i < size; i++)
      //      BaryonField[Vel1Num+dim][i] += RandomForcingField[dim][i]*levelNorm;
      BaryonField[Vel1Num+dim][i] += (RandomForcingField[dim][i]-bulkMomentum[dim])*levelNorm;
 
  return SUCCESS;
 
}
