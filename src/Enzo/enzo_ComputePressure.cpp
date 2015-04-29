// See LICENSE_ENZO file for license and copyright information

/***********************************************************************
/
/  GRID CLASS (COMPUTE THE PRESSURE FIELD AT THE GIVEN TIME)
/
/  written by: Greg Bryan
/  date:       November, 1994
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
 
//----------------------------------------------------------------------
 
int EnzoBlock::ComputePressure(enzo_float time, 
			       enzo_float *pressure,
			       int comoving_coordinates)
{
  /* declarations */
 
  int i, size = 1;
 
  /* Error Check */

  /* Compute the size of the grid. */

  int rank = this->rank();

  for (int dim = 0; dim < rank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  Field field = data()->field();
  
  enzo_float * density         = (enzo_float *) field.values("density");
  enzo_float * total_energy    = (enzo_float *) field.values("total_energy");
  enzo_float * internal_energy = (enzo_float *) field.values("internal_energy");
  enzo_float * velocity_x      = (rank >= 1) ? 
    (enzo_float *) field.values("velocity_x") : NULL;
  enzo_float * velocity_y      = (rank >= 2) ? 
    (enzo_float *) field.values("velocity_y") : NULL;
  enzo_float * velocity_z      = (rank >= 3) ? 
    (enzo_float *) field.values("velocity_z") : NULL;

  /* Loop over the grid, compute the thermal energy, then the pressure,
     the timestep and finally the implied timestep. */
 
  /* special loop for no interpolate. */

  // WARNING: floating point comparison
  if (time == this->time()) {

    if (rank >= 1) {
      for (i = 0; i < size; i++) {
	enzo_float te = total_energy[i];
	enzo_float vx = (rank >= 1) ? velocity_x[i] : 0.0;
	internal_energy[i] = te - 0.5*vx*vx;
      }
    }
    if (rank >= 2) {
      for (i = 0; i < size; i++) {
	enzo_float vy = (rank >= 2) ? velocity_y[i] : 0.0;
	internal_energy[i] -= 0.5*vy*vy;
      }
    }
    if (rank >= 3) {
      for (i = 0; i < size; i++) {
	enzo_float vz = (rank >= 3) ? velocity_z[i] : 0.0;
	internal_energy[i] -= 0.5*vz*vz;
      }
    }

    for (i = 0; i < size; i++) {

      /* gas energy = E - 1/2 v^2. */
 
      enzo_float d  = density[i];

      pressure[i] = (Gamma - 1.0)*d*internal_energy[i];

      if (pressure[i] < pressure_floor)
	pressure[i] = pressure_floor;

    } // end of loop

  } else

    ERROR("EnzoBlock::ComputePressure()",
	  "Accessing OldBaryonField");
	    
  /* general case: */
 
  // for (i = 0; i < size; i++) {
 
  //   enzo_float te = coef * te  BaryonField[TENum][i] +
  // 	              coefold*OldBaryonField[TENum][i];
  //   density       = coef   *   BaryonField[DensNum][i] +
  //                   coefold*OldBaryonField[DensNum][i];
  //   vx     = coef   *   BaryonField[Vel1Num][i] +
  //                   coefold*OldBaryonField[Vel1Num][i];
 
  //   if (GridRank > 1)
  // 	vy   = coef   *   BaryonField[Vel2Num][i] +
  // 	              coefold*OldBaryonField[Vel2Num][i];
  //   if (GridRank > 2)
  // 	vz   = coef   *   BaryonField[Vel3Num][i] +
  // 	              coefold*OldBaryonField[Vel3Num][i];
 
  //   /* gas energy = E - 1/2 v^2. */
 
  //   gas_energy    = total_energy - 0.5*(vx*vx +
  // 					  vy*vy +
  // 					  vz*vz);
 
  //   pressure[i] = (Gamma - 1.0)*density*gas_energy;
 
  //   if (pressure[i] < pressure_floor)
  // 	pressure[i] = pressure_floor;
 
  // }
 
  /* Correct for Gamma from H2. */
 
  if (MultiSpecies > 1) {
 
    enzo_float TemperatureUnits = 1, number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(Gamma-1.0), x, Gamma1, temp;
    enzo_float DensityUnits, LengthUnits, VelocityUnits, TimeUnits;
 
    /* Find Multi-species fields. */

   
    enzo_float * species_De    = (enzo_float *) field.values("species_De");
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

    /* Find the temperature units if we are using comoving coordinates. */
 
    if (comoving_coordinates)
      if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			    &TimeUnits, &VelocityUnits, this->time()) == ENZO_FAIL) {
	fprintf(stderr, "Error in CosmologyGetUnits.\n");
	return ENZO_FAIL;
      }
 
    for (i = 0; i < size; i++) {
 
      number_density =
	0.25*(species_HeI[i] + 
	      species_HeII[i] + 
	      species_HeIII[i]) 
	+    (species_HI[i] +
	      species_HII[i] + 
	      species_De[i]);
 
      nH2 = 0.5*(species_H2I[i] + species_H2II[i]);
 
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
 
    /* To emulate the opacity limit in turbulent star formation 
       simulations */
  
  enzo_float Gamma1 = Gamma;
  if ((ProblemType == 60 || ProblemType == 61))
    for (i=0; i<size; i++) {
      Gamma1 = MIN(Gamma + (log10(density[i])-8.0)*0.3999/2.5, 1.4);
      pressure[i] *= (Gamma1 - 1.0)/(Gamma - 1.0);
    }


  return ENZO_SUCCESS;
}
