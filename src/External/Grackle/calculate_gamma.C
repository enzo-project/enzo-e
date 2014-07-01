/***********************************************************************
/
/ Calculate gamma (ratio of specific heats) field
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle.h"
#include "chemistry_data.h"
#include "code_units.h"
#include "phys_constants.h"

int calculate_temperature(chemistry_data &my_chemistry,
                          code_units &my_units,
                          gr_int grid_rank, gr_int *grid_dimension,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                          gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                          gr_float *H2I_density, gr_float *H2II_density,
                          gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                          gr_float *e_density, gr_float *metal_density,
                          gr_float *temperature);

int calculate_gamma(chemistry_data &my_chemistry,
                    code_units &my_units,
                    gr_int grid_rank, gr_int *grid_dimension,
                    gr_float *density, gr_float *internal_energy,
                    gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                    gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                    gr_float *H2I_density, gr_float *H2II_density,
                    gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                    gr_float *e_density, gr_float *metal_density,
                    gr_float *my_gamma)
{

  if (!my_chemistry.use_grackle)
    return SUCCESS;
 
  gr_int i, size = 1;
  for (int dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];
  
  /* If molecular hydrogen is not being used, just use monotonic.
     (this should not really be called, but provide it just in case). */
 
  for (i = 0; i < size; i++) {
    my_gamma[i] = my_chemistry.Gamma;
  }
 
  if (my_chemistry.primordial_chemistry > 1) {

    /* Compute the temperature first. */
 
    if (calculate_temperature(my_chemistry, my_units,
                              grid_rank, grid_dimension,
                              density, internal_energy,
                              HI_density, HII_density, HM_density,
                              HeI_density, HeII_density, HeIII_density,
                              H2I_density, H2II_density,
                              DI_density, DII_density, HDI_density,
                              e_density, metal_density,
                              my_gamma) == FAIL) {
      fprintf(stderr, "Error in calculate_gamma.\n");
      return FAIL;
    }

    /* Compute Gamma with molecular Hydrogen formula from Omukau \& Nishi
       astro-ph/9811308. */
 
    gr_float x, nH2, number_density, GammaH2Inverse, 
      GammaInverse = 1 / (my_chemistry.Gamma - 1.0);

    for (i = 0; i < size; i++) {
 
      /* Compute relative number abundence of molecular hydrogen. */
 
      number_density =
        0.25 * (HeI_density[i] + HeII_density[i] + HeIII_density[i]) +
        HI_density[i] + HII_density[i] + HM_density[i] +
        e_density[i];
 
      nH2 = 0.5 * (H2I_density[i]  + H2II_density[i]);
 
      /* Only do full computation if there is a reasonable amount of H2.
         The second term in GammaH2Inverse accounts for the vibrational
         degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2 / number_density > 1e-3) {
	x = 6100.0 / my_gamma[i];
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
      }
 
      /* Add in H2. */
 
      my_gamma[i] = 1.0 + (nH2 + number_density) /
        (nH2 * GammaH2Inverse + number_density * GammaInverse);
 
    } // end: loop over i
 
  } // end: if (my_chemistry.primordial_chemistry > 1)

  return SUCCESS;
}
