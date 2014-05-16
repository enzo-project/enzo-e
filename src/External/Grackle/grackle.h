/***********************************************************************
/
/ Grackle function prototypes
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_H__
#define __GRACKLE_H__

#include "grackle_types.h"
#include "code_units.h"
#include "chemistry_data.h"

chemistry_data set_default_chemistry_parameters();

int initialize_chemistry_data(chemistry_data &my_chemistry,
                              code_units &my_units, gr_float a_value);

int initialize_UVbackground_data(chemistry_data &my_chemistry);
int update_UVbackground_rates(chemistry_data &my_chemistry,
			      code_units &my_units, gr_float a_value);

int solve_chemistry(chemistry_data &my_chemistry,
		    code_units &my_units,
		    gr_float a_value, gr_float dt_value,
		    gr_int grid_rank, gr_int *grid_dimension,
		    gr_int *grid_start, gr_int *grid_end,
		    gr_float *density, gr_float *internal_energy,
		    gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
		    gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
		    gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
		    gr_float *H2I_density, gr_float *H2II_density,
		    gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
		    gr_float *e_density, gr_float *metal_density);

int calculate_cooling_time(chemistry_data &my_chemistry,
			   code_units &my_units, gr_float a_value,
			   gr_int grid_rank, gr_int *grid_dimension,
			   gr_int *grid_start, gr_int *grid_end,
			   gr_float *density, gr_float *internal_energy,
			   gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
			   gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
			   gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
			   gr_float *H2I_density, gr_float *H2II_density,
			   gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
			   gr_float *e_density, gr_float *metal_density,
			   gr_float *cooling_time);

int calculate_gamma(chemistry_data &my_chemistry,
                    code_units &my_units,
                    gr_int grid_rank, gr_int *grid_dimension,
                    gr_float *density, gr_float *internal_energy,
                    gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                    gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                    gr_float *H2I_density, gr_float *H2II_density,
                    gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                    gr_float *e_density, gr_float *metal_density,
                    gr_float *my_gamma);

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
                          gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                          gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                          gr_float *H2I_density, gr_float *H2II_density,
                          gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                          gr_float *e_density, gr_float *metal_density,
                          gr_float *temperature);

/* Tabular-only functions. */

int solve_chemistry(chemistry_data &my_chemistry,
                    code_units &my_units,
                    gr_float a_value, gr_float dt_value,
                    gr_int grid_rank, gr_int *grid_dimension,
                    gr_int *grid_start, gr_int *grid_end,
                    gr_float *density, gr_float *internal_energy,
                    gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                    gr_float *metal_density);

int calculate_cooling_time(chemistry_data &my_chemistry,
                           code_units &my_units, gr_float a_value,
                           gr_int grid_rank, gr_int *grid_dimension,
                           gr_int *grid_start, gr_int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                           gr_float *metal_density,
                           gr_float *cooling_time);

int calculate_pressure(chemistry_data &my_chemistry,
                       code_units &my_units,
                       gr_int grid_rank, gr_int *grid_dimension,
                       gr_float *density, gr_float *internal_energy,
                       gr_float *pressure);

int calculate_temperature(chemistry_data &my_chemistry,
                          code_units &my_units,
                          gr_int grid_rank, gr_int *grid_dimension,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *metal_density,
                          gr_float *temperature);

#endif
