/***********************************************************************
/
/ Calculate temperature field (tabulated cooling function)
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
 
/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
  
/* function prototypes */ 
 
int calculate_temperature(chemistry_data &my_chemistry,
                          code_units &my_units,
                          gr_int grid_rank, gr_int *grid_dimension,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *metal_density,
                          gr_float *temperature)
{

  if (!my_chemistry.use_grackle)
    return SUCCESS;

  if (my_chemistry.primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }
 
  /* Compute the size of the fields. */
 
  gr_int i, size = 1;
  for (int dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  /* Calculate temperature units. */

  gr_float temperature_units =  mh * POW(my_units.velocity_units, 2) / kboltz;

  gr_float number_density, tiny_number = 1.-20;
  gr_float inv_metal_mol = 1.0 / MU_METAL;

  gr_float munew, muold;
  gr_int ti, ti_max, index;
  ti_max = 20;

  gr_float logtem0 = log(my_chemistry.TemperatureStart);
  gr_float logtem9 = log(my_chemistry.TemperatureEnd);
  gr_float dlogtem = (log(my_chemistry.TemperatureEnd) - 
                      log(my_chemistry.TemperatureStart)) / 
                      (my_chemistry.NumberOfTemperatureBins - 1);
  gr_float logtem, t1, t2, tdef;

  /* Compute temperature with mu calculated directly. */
 
  for (i = 0; i < size; i++) {

    munew = 1.0;

    for (ti = 0; ti < ti_max; ti++) {

      muold = munew;
      temperature[i] = std::max((my_chemistry.Gamma - 1.) * 
                           internal_energy[i] *
                           munew * temperature_units,
                           my_chemistry.TemperatureStart);
      logtem = log(temperature[i]);
      logtem = std::max(logtem, logtem0);
      logtem = std::min(logtem, logtem9);

      index = std::min((const int) my_chemistry.NumberOfTemperatureBins - 2,
                  std::max(0, (int) ((logtem-logtem0)/dlogtem)));
      t1 = (logtem0 + (index)     * dlogtem);
      t2 = (logtem0 + (index + 1) * dlogtem);
      tdef = (logtem - t1) / (t2 - t1);
      munew = my_chemistry.mu[index] + tdef
        * (my_chemistry.mu[index+1] - my_chemistry.mu[index]);

      temperature[i] = temperature[i] * munew / muold;

      if (fabs((munew/muold) - 1.) <= 1.e-2) {
        muold = munew;

        // Add metal species to mean molecular weight
          
        munew = density[i] / (density[i] / munew +
                              metal_density[i] / MU_METAL);
        temperature[i] = temperature[i] * munew / muold;
        break;

      }

    }

    if (ti >= ti_max) {
      fprintf(stderr, "Warning: mean molecular weight failed to converge!\n");
    }

  }

  return SUCCESS;
}
