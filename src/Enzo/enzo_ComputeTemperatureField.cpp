// See LICENSE_ENZO file for license and copyright information

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

#include "cello.hpp"

#include "enzo.hpp"
 
/* Set the mean molecular mass. */
 
#define DEFAULT_MU 0.6
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
//---------------------------------------------------------------------- 
 
int EnzoBlock::ComputeTemperatureField
(
 enzo_float *temperature,
 int comoving_coordinates
 )
{
 
  int result;

  Field field = data()->field();

  /* If Gadget equilibrium cooling is on, call the appropriate routine,
     then exit - don't use the rest of the routine. */

//   if(GadgetEquilibriumCooling){
//     if(DualEnergyFormalism)
//       result = this->GadgetComputeTemperatureDEF(Time(), temperature);
//     else
//       result = this->GadgetComputeTemperature(Time(),temperature);

//     if(result == ENZO_FAIL) {
//       fprintf(stderr, "Error in grid->ComputePressure: Gadget.\n");
//       return ENZO_FAIL;
//     }
//     return ENZO_SUCCESS;
//   }

  /* Compute the pressure first. */
 
  if (DualEnergyFormalism)
    result = ComputePressureDualEnergyFormalism(time(), 
						temperature,
						comoving_coordinates);
  else
    result = ComputePressure(time(), temperature, comoving_coordinates);
 
  if (result == ENZO_FAIL) {
    fprintf(stderr, "Error in grid->ComputePressure.\n");
    return ENZO_FAIL;
  }
 
  /* Compute the size of the fields. */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  enzo_float * density = (enzo_float *) field.values("density");

  // @@@ WHY PROBLEM-DEPENDENT? jb @@@
  if (ProblemType == 60 || ProblemType == 61) { //AK
    for (i = 0; i < size; i++) {
      if (density[i] <= 0.0)
	temperature[i] = 1.0;
      else
	temperature[i] /=  density[i];
    }
    return ENZO_SUCCESS;
  }
 
  enzo_float TemperatureUnits = 1, number_density;
  enzo_float DensityUnits, LengthUnits, VelocityUnits, TimeUnits;
 
  /* Find the temperature units if we are using comoving coordinates. */
 
  if (comoving_coordinates)
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, time()) == ENZO_FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return ENZO_FAIL;
    }

  /* For Sedov Explosion compute temperature without floor */

  // @@@ WHY PROBLEM-DEPENDENT? jb @@@
  enzo_float mol_weight = DEFAULT_MU, min_temperature = 1.0;
  if (ProblemType == 7) {//AK for Sedov explosion test
    mol_weight = 1.0;
    min_temperature = temperature_floor;
  }

  if (MultiSpecies == FALSE)
 
    /* If the multi-species flag is not set,
       Compute temperature T = p/d and assume mu = DEFAULT_MU. */
 
    for (i = 0; i < size; i++)
      temperature[i] = MAX((TemperatureUnits*temperature[i]*mol_weight
			    /MAX(density[i], 
				 (enzo_float)(density_floor))),
			 min_temperature);
  else {
 
    /* Find Multi-species fields. */

    enzo_float * species_HI    = (enzo_float *) field.values("species_HI");
    enzo_float * species_HII   = (enzo_float *) field.values("species_HII");
    enzo_float * species_HeI   = (enzo_float *) field.values("species_HeI");
    enzo_float * species_HeII  = (enzo_float *) field.values("species_HeII");
    enzo_float * species_HeIII = (enzo_float *) field.values("species_HeIII");
    enzo_float * species_HM    = (enzo_float *) field.values("species_HM");
    enzo_float * species_H2I   = (enzo_float *) field.values("species_H2I");
    enzo_float * species_H2II  = (enzo_float *) field.values("species_H2II");
    enzo_float * species_DI    = (enzo_float *) field.values("species_DI");
    enzo_float * species_DII   = (enzo_float *) field.values("species_DII");
    enzo_float * species_HDI   = (enzo_float *) field.values("species_HDI");
 
    /* Compute temperature with mu calculated directly. */
 
    for (i = 0; i < size; i++) {
 
      number_density =
	0.25*(species_HeI[i]  + 
	      species_HeII[i] +
	      species_HeIII[i] ) +
	species_HI[i] + species_HII[i] + density[i];
 
      /* Add in H2. */
 
      if (MultiSpecies > 1)
	number_density += species_HM[i]   +
	  0.5*(species_H2I[i] +
	       species_H2II[i]);
 
      /* Ignore deuterium. */
 
      temperature[i] *= TemperatureUnits/
	MAX(number_density,  number_density_floor);
      temperature[i] = MAX(temperature[i], 
			   MINIMUM_TEMPERATURE);
    }
  }
 
  return ENZO_SUCCESS;
}
