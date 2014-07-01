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
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */ 

int calculate_pressure(chemistry_data &my_chemistry,
                       code_units &my_units,
                       gr_int grid_rank, gr_int *grid_dimension,
                       gr_float *density, gr_float *internal_energy,
                       gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                       gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                       gr_float *H2I_density, gr_float *H2II_density,
                       gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                       gr_float *e_density, gr_float *metal_density,
                       gr_float *pressure);

int calculate_temperature(chemistry_data &my_chemistry,
                          code_units &my_units,
                          gr_int grid_rank, gr_int *grid_dimension,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *metal_density,
                          gr_float *temperature);
 
int calculate_temperature(chemistry_data &my_chemistry,
                          code_units &my_units,
                          gr_int grid_rank, gr_int *grid_dimension,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                          gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                          gr_float *H2I_density, gr_float *H2II_density,
                          gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                          gr_float *e_density, gr_float *metal_density,
                          gr_float *temperature)
{

  if (!my_chemistry.use_grackle)
    return SUCCESS;

  /* Compute the pressure first. */
 
  if (calculate_pressure(my_chemistry, my_units,
                         grid_rank, grid_dimension,
                         density, internal_energy,
                         HI_density, HII_density, HM_density,
                         HeI_density, HeII_density, HeIII_density,
                         H2I_density, H2II_density,
                         DI_density, DII_density, HDI_density,
                         e_density, metal_density,
                         temperature) == FAIL) {
    fprintf(stderr, "Error in calculate_pressure.\n");
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
  
  if (my_chemistry.primordial_chemistry == 0) {
    if (calculate_temperature(my_chemistry,
                              my_units,
                              grid_rank, grid_dimension,
                              density, internal_energy,
                              metal_density,
                              temperature) == FAIL) {
      fprintf(stderr, "Error in calculcate_temperature_table.\n");
      return FAIL;
    }
    return SUCCESS;
  }

 /* Compute temperature with mu calculated directly. */
 
  for (i = 0; i < size; i++) {
 
    if (my_chemistry.primordial_chemistry > 0) {
      number_density =
        0.25 * (HeI_density[i] + HeII_density[i] +  HeIII_density[i]) +
        HI_density[i] + HII_density[i] + e_density[i];
    }

    /* Add in H2. */
 
    if (my_chemistry.primordial_chemistry > 1) {
      number_density += HM_density[i] + 
        0.5 * (H2I_density[i] + H2II_density[i]);
    }

    if (metal_density != NULL) {
      number_density += metal_density[i] * inv_metal_mol;
    }
 
    /* Ignore deuterium. */
 
    temperature[i] *= temperature_units / std::max(number_density, tiny_number);
    temperature[i] = std::max(temperature[i], MINIMUM_TEMPERATURE);
  }
 
  return SUCCESS;
}
