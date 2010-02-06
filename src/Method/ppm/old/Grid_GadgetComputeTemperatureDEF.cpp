/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * SYNOPSIS:
 *
 *    
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
/***********************************************************************
/
/  GRID CLASS (COMPUTE THE TEMPERATURE AT THE GIVEN TIME USING
/              GADGET EQUILIBRIUM COOLING - DUAL ENERGY
/
/  written by: Brian O'Shea
/  date:       February 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
/  NOTE:  This routine is pretty blatantly borrowed from 
/         Grid_ComputePressureDualEnergyFormalism - if there are bugs 
/         here, check there as well.
/
************************************************************************/

// Compute the pressure at the requested time.  The pressure here is
//   just the ideal-gas equation-of-state (dual energy version).

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int grid::GadgetComputeTemperatureDEF(FLOAT time, float *temperature)
{

  /* declarations */

  float density, gas_energy;
  int i, size = 1;

  float *ne_guess,guess;
  guess=0.01;
  ne_guess = &guess;

  /* Error Check */

  if (time < OldTime || time > Time) {
    fprintf(stderr, "requested time is outside available range.\n");
    return FAIL;
  }

  /* Compute interpolation coefficients. */

  float coef, coefold;
  if (Time - OldTime > 0)
    coef    = (time - OldTime)/(Time - OldTime);
  else
    coef    = 1;

  coefold = 1 - coef;

  /* Compute the size of the grid. */

  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  // get physical units
  float DensityUnits, LengthUnits, VelocityUnits, TimeUnits,TemperatureUnits=1;

  if (ComovingCoordinates)
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }

  /* Loop over the grid, compute the thermal energy, then the temperature,
     the timestep and finally the implied timestep. */

  /* special loop for no interpolate. */

  if (time == Time)

    for (i = 0; i < size; i++) {
      gas_energy = BaryonField[GENum][i];
      density = BaryonField[DensNum][i];
      
      gas_energy *= (VelocityUnits*VelocityUnits);

      density *= DensityUnits;

      temperature[i] = Gadgetconvert_u_to_temp(gas_energy, density, ne_guess);
      
      if (temperature[i] < 1.0)
	temperature[i] = 1.0;
      
    }

  else

    /* general case: */

    for (i = 0; i < size; i++) {

      gas_energy    = coef   *   BaryonField[GENum][i] + 
	              coefold*OldBaryonField[GENum][i];
      density       = coef   *   BaryonField[DensNum][i] + 
                      coefold*OldBaryonField[DensNum][i];


      gas_energy *= (VelocityUnits*VelocityUnits);

      density *= DensityUnits;

      temperature[i] = Gadgetconvert_u_to_temp(gas_energy, density, ne_guess);

      if (temperature[i] < 1.0)
	temperature[i] = 1.0;



    }

  return SUCCESS;
}
