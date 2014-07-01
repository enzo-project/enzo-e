/***********************************************************************
/
/ Calculate pressure field (tabulated cooling function)
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
                       gr_float *pressure)
{

  if (!my_chemistry.use_grackle)
    return SUCCESS;

  if (my_chemistry.primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }

  gr_float tiny_number = 1.e-20;
  gr_int i, size = 1;
  for (int dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  for (i = 0; i < size; i++) {
 
    pressure[i] = (my_chemistry.Gamma - 1.0) * density[i] * internal_energy[i];
 
    if (pressure[i] < tiny_number)
      pressure[i] = tiny_number;
  }  
 
  return SUCCESS;
}
