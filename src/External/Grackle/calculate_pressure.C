/***********************************************************************
/
/ Calculate pressure field
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

int calculate_pressure(chemistry_data &my_chemistry,
                       code_units &my_units,
                       gr_int grid_rank, gr_int *grid_dimension,
                       gr_float *density, gr_float *internal_energy,
                       gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                       gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                       gr_float *H2I_density, gr_float *H2II_density,
                       gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                       gr_float *e_density, gr_float *metal_density,
                       gr_float *pressure)
{

  if (!my_chemistry.use_grackle)
    return SUCCESS;

  gr_float tiny_number = 1.e-20;
  gr_int i, size = 1;
  for (int dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  for (i = 0; i < size; i++) {
 
    pressure[i] = (my_chemistry.Gamma - 1.0) * density[i] * internal_energy[i];
 
    if (pressure[i] < tiny_number)
      pressure[i] = tiny_number;
  }
  
  /* Correct for Gamma from H2. */
 
  if (my_chemistry.primordial_chemistry > 1) {
 
    /* Calculate temperature units. */

    gr_float temperature_units =  mh * POW(my_units.velocity_units, 2) / kboltz;

    gr_float number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(my_chemistry.Gamma-1.0), x, Gamma1, temp;
  
    for (i = 0; i < size; i++) {
 
      number_density =
        0.25 * (HeI_density[i] + HeII_density[i] + HeIII_density[i]) +
        HI_density[i] + HII_density[i] + HM_density[i] +
        e_density[i];
 
      nH2 = 0.5 * (H2I_density[i] + H2II_density[i]);
 
      /* First, approximate temperature. */
 
      if (number_density == 0)
	number_density = tiny_number;
      temp = std::max(temperature_units * pressure[i] / (number_density + nH2), 1.0);
 
      /* Only do full computation if there is a reasonable amount of H2.
	 The second term in GammaH2Inverse accounts for the vibrational
	 degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2 / number_density > 1e-3) {
        x = 6100.0 / temp;
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
      }
 
      Gamma1 = 1.0 + (nH2 + number_density) /
	             (nH2 * GammaH2Inverse + number_density * GammaInverse);
	
      /* Correct pressure with improved Gamma. */
 
      pressure[i] *= (Gamma1 - 1.0) / (my_chemistry.Gamma - 1.0);
 
    } // end: loop over i
 
  } // end: if (my_chemistry.primordial_chemistry > 1)
 
  return SUCCESS;
}
