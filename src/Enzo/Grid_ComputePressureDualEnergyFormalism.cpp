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
 
  int i, size = 1;
 
  /* Compute the size of the grid. */
 
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  Field field = block()->field();

  enzo_float * density         = (enzo_float *) field.values("density");
  enzo_float * internal_energy = (enzo_float *) field.values("internal_energy");

  /* Loop over the grid, compute the thermal energy, then the pressure,
     the timestep and finally the implied timestep. */
 
  /* special loop for no interpolate. */
 
  if (time == this->time()) {
 
    for (i = 0; i < size; i++) {
      pressure[i] = (Gamma - 1.0) * density[i] *
                                    internal_energy[i];
 
      if (pressure[i] < pressure_floor)
	pressure[i] = pressure_floor;
    }
 
  } else {
 
    /* general case: */

    ERROR("EnzoBlock::ComputePressure()",
	    "Accessing OldBaryonField");

    // for (i = 0; i < size; i++) {
 
    //   gas_energy    = coef   *   internal_energy[i] +
    // 	coefold*OldBaryonField[GENum][i];
    //   density       = coef   *   density[i] +
    // 	coefold*OldBaryonField[DensNum][i][i];
 
    //   pressure[i] = (Gamma - 1.0)*density*gas_energy;
 
    //   if (pressure[i] < pressure_floor)
    // 	pressure[i] = pressure_floor;
 
    // }
  }

  /* Correct for Gamma from H2. */
 
  if (MultiSpecies > 1) {
 
    enzo_float TemperatureUnits = 1, number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(Gamma-1.0), x, Gamma1, temp;
    enzo_float DensityUnits, LengthUnits, VelocityUnits, TimeUnits;
 
    enzo_float * species_De    = (enzo_float *) field.values("species_De");
    enzo_float * species_HI    = (enzo_float *) field.values("species_HI");
    enzo_float * species_HII   = (enzo_float *) field.values("species_HII");
    enzo_float * species_HeI   = (enzo_float *) field.values("species_HeI");
    enzo_float * species_HeII  = (enzo_float *) field.values("species_HeII");
    enzo_float * species_HeIII = (enzo_float *) field.values("species_HeIII");
    // enzo_float * species_HM    = (enzo_float *) field.values("species_HM");
    enzo_float * species_H2I   = (enzo_float *) field.values("species_H2I");
    enzo_float * species_H2II  = (enzo_float *) field.values("species_H2II");
    // enzo_float * species_DI    = (enzo_float *) field.values("species_DI");
    // enzo_float * species_DII   = (enzo_float *) field.values("species_DII");
    // enzo_float * species_HDI   = (enzo_float *) field.values("species_HDI");
 
    /* Find the temperature units if we are using comoving coordinates. */
 
    if (ComovingCoordinates)
      if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			    &TimeUnits, &VelocityUnits, this->time()) == ENZO_FAIL) {
	fprintf(stderr, "Error in CosmologyGetUnits.\n");
	return ENZO_FAIL;
      }
 
    for (i = 0; i < size; i++) {
 
      number_density =
	  0.25 * (species_HeI[i]  + 
		  species_HeII[i] +
		  species_HeIII[i] )
	+        (species_HI[i]   + 
		  species_HII[i]  +
		  species_De[i]);
 
      nH2 = 0.5*(species_H2I[i]  + species_H2II[i]);
 
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
