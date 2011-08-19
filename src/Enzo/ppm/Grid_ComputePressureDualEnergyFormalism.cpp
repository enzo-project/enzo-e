// See LICENSE_ENZO file for license and copyright information

/// @file      Grid_ComputePressureDualEnergyFormalism.cpp
/// @author    Greg Bryan
/// @date      November, 1994
/// @brief     (COMPUTE THE PRESSURE FIELD AT THE GIVEN TIME) - DUAL ENERGY
///
/// Compute the pressure at the requested time.  The pressure here is
/// just the ideal-gas equation-of-state (dual energy version).

#include "cello.hpp"

#include "enzo.hpp"
 
//----------------------------------------------------------------------
 
int EnzoBlock::ComputePressureDualEnergyFormalism
(enzo_float time, enzo_float *pressure)
{
 
  /* declarations */
 
  enzo_float density, gas_energy;
  int i, size = 1;
 
  /* Error Check */
 
  if (time < OldTime || time > Time()) {
    fprintf(stderr, "requested time is outside available range.\n");
    return ENZO_FAIL;
  }
 
  /* Compute interpolation coefficients. */
 
  enzo_float coef, coefold;
  if (Time() - OldTime > 0)
    coef    = (time - OldTime)/(Time() - OldTime);
  else
    coef    = 1;
 
  coefold = 1 - coef;
 
  /* Compute the size of the grid. */
 
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == ENZO_FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return ENZO_FAIL;
  }
 
  /* Loop over the grid, compute the thermal energy, then the pressure,
     the timestep and finally the implied timestep. */
 
  /* special loop for no interpolate. */
 
  if (time == Time())
 
    for (i = 0; i < size; i++) {
      pressure[i] = (Gamma - 1.0) * BaryonField[DensNum][i] *
                                    BaryonField[GENum][i];
 
      if (pressure[i] < pressure_floor)
	pressure[i] = pressure_floor;
 
    }
 
  else
 
    /* general case: */
 
    for (i = 0; i < size; i++) {
 
      gas_energy    = coef   *   BaryonField[GENum][i] +
	              coefold*OldBaryonField[GENum][i];
      density       = coef   *   BaryonField[DensNum][i] +
                      coefold*OldBaryonField[DensNum][i];
 
      pressure[i] = (Gamma - 1.0)*density*gas_energy;
 
      if (pressure[i] < pressure_floor)
	pressure[i] = pressure_floor;
 
    }
 
  /* Correct for Gamma from H2. */
 
  if (MultiSpecies > 1) {
 
    enzo_float TemperatureUnits = 1, number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(Gamma-1.0), x, Gamma1, temp;
    enzo_float DensityUnits, LengthUnits, VelocityUnits, TimeUnits;
 
    /* Find Multi-species fields. */
 
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
        H2IINum, DINum, DIINum, HDINum;
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == ENZO_FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return ENZO_FAIL;
    }
 
    /* Find the temperature units if we are using comoving coordinates. */
 
    if (ComovingCoordinates)
      if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			    &TimeUnits, &VelocityUnits, Time()) == ENZO_FAIL) {
	fprintf(stderr, "Error in CosmologyGetUnits.\n");
	return ENZO_FAIL;
      }
 
    for (i = 0; i < size; i++) {
 
      number_density =
	  0.25*(BaryonField[HeINum][i]  + BaryonField[HeIINum][i] +
		BaryonField[HeIIINum][i]                        ) +
	        BaryonField[HINum][i]   + BaryonField[HIINum][i]  +
                BaryonField[DeNum][i];
 
      nH2 = 0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);
 
      /* First, approximate temperature. */
 
      if (number_density == 0)
	number_density = number_density_floor;
      temp = MAX(TemperatureUnits*pressure[i]/(number_density + nH2), 
		 (enzo_float)(1.0));
 
      /* Only do full computation if there is a reasonable amount of H2.
	 The second term in GammaH2Inverse accounts for the vibrational
	 degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2/number_density > 1e-3) {
	x = temp/6100.0;
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/pow(exp(x)-1.0,2));
      }
 
      Gamma1 = 1.0 + (nH2 + number_density) /
	             (nH2*GammaH2Inverse + number_density * GammaInverse);
	
      /* Correct pressure with improved Gamma. */
 
      pressure[i] *= (Gamma1 - 1.0)/(Gamma - 1.0);
 
    } // end: loop over i
 
  } // end: if (MultiSpecies > 1)
 
  return ENZO_SUCCESS;
}
