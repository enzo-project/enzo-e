// See LICENSE_ENZO file for license and copyright information

/***********************************************************************
/
/  GRID CLASS (COMPUTE THE GAMMA FIELD)
/
/  written by: Greg Bryan
/  date:       May, 1999
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the ratio of specific heats.

#include "cello.hpp"

#include "enzo.hpp"
 
/* function prototypes */

int EnzoBlock::ComputeGammaField(enzo_float *GammaField)
{
 
  /* Compute the size of the fields. */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  if (MultiSpecies < 2)
 
    /* If molecular hydrogen is not being used, just use monotonic.
       (this should not really be called, but provide it just in case). */
 
    for (i = 0; i < size; i++)
      GammaField[i] = Gamma;
 
  else {
 

    Field field = block()->field();
   
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

    /* Compute the temperature. */

    ComputeTemperatureField(GammaField);
 
    /* Compute Gamma with molecular Hydrogen formula from Omukau \& Nishi
       astro-ph/9811308. */
 
    enzo_float x, nH2, number_density, GammaH2Inverse, GammaInverse = 1/(Gamma-1.0);
    for (i = 0; i < size; i++) {
 
      /* Compute relative number abundence of molecular hydrogen. */
 
      number_density =
	0.25*(species_HeI[i]  + species_HeII[i] +
	      species_HeIII[i]                        ) +
              species_HI[i]   + species_HII[i]  +
              species_De[i];
 
      nH2 = 0.5*(species_H2I[i]  + species_H2II[i]);
 
      /* Only do full computation if there is a reasonable amount of H2.
         The second term in GammaH2Inverse accounts for the vibrational
         degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2/number_density > 1e-3) {
	x = GammaField[i]/6100.0;
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/pow(exp(x)-1.0,2));
      }
 
      /* Add in H2. */
 
      GammaField[i] = 1.0 + (nH2 + number_density) /
                      (nH2*GammaH2Inverse + number_density * GammaInverse);
 
    }
 
  } // end: if (MultiSpecies < 2)
 
  return ENZO_SUCCESS;
}
