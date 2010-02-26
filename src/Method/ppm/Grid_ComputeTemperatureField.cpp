/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
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
/  GRID CLASS (COMPUTE THE TEMPERATURE FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the pressure at the requested time.  The pressure here is
//   just the ideal-gas equation-of-state.

#include "cello_hydro.h"
 
/* Set the mean molecular mass. */
 
#define DEFAULT_MU 0.6
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
 
 
int ComputeTemperatureField(float *temperature)
{
 
  int DensNum, result;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
 
  /* If Gadget equilibrium cooling is on, call the appropriate routine,
     then exit - don't use the rest of the routine. */

//   if(GadgetEquilibriumCooling){
//     if(DualEnergyFormalism)
//       result = this->GadgetComputeTemperatureDEF(Time, temperature);
//     else
//       result = this->GadgetComputeTemperature(Time,temperature);

//     if(result == FAIL) {
//       fprintf(stderr, "Error in grid->ComputePressure: Gadget.\n");
//       return FAIL;
//     }
//     return SUCCESS;
//   }

  /* Compute the pressure first. */
 
  if (DualEnergyFormalism)
    result = ComputePressureDualEnergyFormalism(Time, temperature);
  else
    result = ComputePressure(Time, temperature);
 
  if (result == FAIL) {
    fprintf(stderr, "Error in grid->ComputePressure.\n");
    return FAIL;
  }
 
  /* Compute the size of the fields. */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find density.\n");
    return FAIL;
  }
 

  // @@@ WHY PROBLEM-DEPENDENT? jb @@@
  if (ProblemType == 60 || ProblemType == 61) { //AK
    for (i = 0; i < size; i++) {
      if (BaryonField[DensNum][i] <= 0.0)
	temperature[i] = 1.0;
      else
	temperature[i] /=  BaryonField[DensNum][i];
    }
    return SUCCESS;
  }
 
  float TemperatureUnits = 1, number_density;
  float DensityUnits, LengthUnits, VelocityUnits, TimeUnits;
 
  /* Find the temperature units if we are using comoving coordinates. */
 
  if (ComovingCoordinates)
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }

  /* For Sedov Explosion compute temperature without floor */

  // @@@ WHY PROBLEM-DEPENDENT? jb @@@
  float mol_weight = DEFAULT_MU, min_temperature = 1.0;
  if (ProblemType == 7) {//AK for Sedov explosion test
    mol_weight = 1.0;
    min_temperature = temperature_floor;
  }


  if (MultiSpecies == FALSE)
 
    /* If the multi-species flag is not set,
       Compute temperature T = p/d and assume mu = DEFAULT_MU. */
 
    for (i = 0; i < size; i++)
      temperature[i] = max((TemperatureUnits*temperature[i]*mol_weight
		         /max(BaryonField[DensNum][i], density_floor)),
			 min_temperature);
  else {
 
    /* Find Multi-species fields. */
 
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return FAIL;
    }
 
    /* Compute temperature with mu calculated directly. */
 
    for (i = 0; i < size; i++) {
 
      number_density =
	0.25*(BaryonField[HeINum][i]  + BaryonField[HeIINum][i] +
	      BaryonField[HeIIINum][i]                        ) +
              BaryonField[HINum][i]   + BaryonField[HIINum][i]  +
              BaryonField[DeNum][i];
 
      /* Add in H2. */
 
      if (MultiSpecies > 1)
	number_density += BaryonField[HMNum][i]   +
	  0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);
 
      /* Ignore deuterium. */
 
      temperature[i] *= TemperatureUnits/max(number_density, number_density_floor);
      temperature[i] = max(temperature[i], MINIMUM_TEMPERATURE);
    }
  }
 
  return SUCCESS;
}
