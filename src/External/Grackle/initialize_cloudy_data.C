/***********************************************************************
/
/ Initialize Cloudy cooling data
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <math.h>
#include "hdf5.h"
#include "grackle_macros.h"
#include "grackle_types.h"
#include "chemistry_data.h"
#include "code_units.h"

#define SMALL_LOG_VALUE -99.0
#define CLOUDY_MAX_DIMENSION 3

/**************************** Functions Prototypes ******************************/

// Initialize Cloudy cooling data
int initialize_cloudy_data(chemistry_data &my_chemistry,
                           cloudy_data &my_cloudy, char *group_name,
                           code_units &my_units, gr_float a_value,
                           gr_int read_data)
{

  gr_int q, w;
  double *temp_data;
  long long temp_int;
  long long *temp_int_arr;
  char parameter_name[MAX_LINE_LENGTH];
  gr_int debug = 0;

  // Initialize things needed even if cloudy cooling is not used.

  my_cloudy.grid_parameters = new gr_float*[CLOUDY_MAX_DIMENSION];
  my_cloudy.grid_dimension = new gr_int[CLOUDY_MAX_DIMENSION];
  for (q = 0;q < CLOUDY_MAX_DIMENSION;q++) {
    my_cloudy.grid_dimension[q] = 0;
  }

  // Zero arrays if cloudy cooling not used.

  if (!read_data) {
    my_cloudy.grid_rank = 0;
    return SUCCESS;
  }

  fprintf(stderr,"Initializing Cloudy cooling: %s.\n", group_name);
  fprintf(stderr,"cloudy_table_file: %s.\n",my_chemistry.grackle_data_file);

  /* Get conversion units. */

  gr_float co_length_units, co_density_units;
  if (my_units.comoving_coordinates == TRUE) {
    co_length_units = my_units.length_units;
    co_density_units = my_units.density_units;
  }
  else {
    co_length_units = my_units.length_units *
      a_value * my_units.a_units;
    co_density_units = my_units.density_units /
      POW(a_value * my_units.a_units, 3);
  }

  double tbase1 = my_units.time_units;
  double xbase1 = co_length_units/(a_value * my_units.a_units);
  double dbase1 = co_density_units * POW(a_value * my_units.a_units, 3);
  double mh = 1.67e-24;
  double CoolUnit = (POW(my_units.a_units,5) * POW(xbase1,2) * POW(mh,2)) /
                    (POW(tbase1,3) * dbase1);

  // Read cooling data in from hdf5 file.

  hid_t       file_id, dset_id, attr_id; 
  herr_t      status;
  herr_t      h5_error = -1;

  file_id = H5Fopen(my_chemistry.grackle_data_file, 
                    H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open cooling dataset and get grid dimensions.

  sprintf(parameter_name, "/CoolingRates/%s/Cooling", group_name);
  dset_id =  H5Dopen(file_id, parameter_name);
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open Cooling in %s.\n",my_chemistry.grackle_data_file);
    return FAIL;
  }

  // Grid rank.
  attr_id = H5Aopen_name(dset_id, "Rank");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open Rank attribute in Cooling dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8, &temp_int);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to read Rank attribute in Cooling dataset.\n");
    return FAIL;
  }
  my_cloudy.grid_rank = (int) temp_int;
  fprintf(stderr,"Cloudy cooling grid rank: %"ISYM".\n",my_cloudy.grid_rank);
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Rank attribute in Cooling dataset.\n");
    return FAIL;
  }

  // Grid dimension.
  temp_int_arr = new long long[my_cloudy.grid_rank];
  attr_id = H5Aopen_name(dset_id, "Dimension");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8,temp_int_arr);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to read Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  fprintf(stderr,"Cloudy cooling grid dimensions:");
  for (q = 0;q < my_cloudy.grid_rank;q++) {
    my_cloudy.grid_dimension[q] = (int) temp_int_arr[q];
    fprintf(stderr," %"ISYM,my_cloudy.grid_dimension[q]);
  }
  fprintf(stderr,".\n");
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  delete [] temp_int_arr;

  // Grid parameters.
  for (q = 0;q < my_cloudy.grid_rank;q++) {

    if (q < my_cloudy.grid_rank - 1) {
      sprintf(parameter_name,"Parameter%"ISYM,(q+1));
    }
    else {
      sprintf(parameter_name,"Temperature");
    }

    temp_data = new double[my_cloudy.grid_dimension[q]];

    attr_id = H5Aopen_name(dset_id, parameter_name);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open %s attribute in Cooling dataset.\n", 
              parameter_name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_R8, temp_data);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to read %s attribute in Cooling dataset.\n",
              parameter_name);
      return FAIL;
    }

    my_cloudy.grid_parameters[q] = 
      new gr_float[my_cloudy.grid_dimension[q]];
    for (w = 0;w < my_cloudy.grid_dimension[q];w++) {
      if (q < my_cloudy.grid_rank - 1) {
	my_cloudy.grid_parameters[q][w] = (float) temp_data[w];
      }
      else {
	// convert temeperature to log
	my_cloudy.grid_parameters[q][w] = (float) log10(temp_data[w]);
      }

    }
    fprintf(stderr,"%s: %"GSYM" to %"GSYM" (%"ISYM" steps).\n",parameter_name,
            my_cloudy.grid_parameters[q][0],
            my_cloudy.grid_parameters[q][my_cloudy.grid_dimension[q]-1],
            my_cloudy.grid_dimension[q]);
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close %s attribute in Cooling dataset.\n",
              parameter_name);
      return FAIL;
    }
    delete [] temp_data;

  }

  // Read Cooling data.
  my_cloudy.data_size = 1;
  for (q = 0;q < my_cloudy.grid_rank;q++) {
    my_cloudy.data_size *= my_cloudy.grid_dimension[q];
  }
  temp_data = new double[my_cloudy.data_size];

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
  fprintf(stderr,"Reading Cloudy Cooling dataset.\n");
  if (status == h5_error) {
    fprintf(stderr,"Failed to read Cooling dataset.\n");
    return FAIL;
  }

  my_cloudy.cooling_data = new gr_float[my_cloudy.data_size];
  for (q = 0;q < my_cloudy.data_size;q++) {
    my_cloudy.cooling_data[q] = temp_data[q] > 0 ? (float) log10(temp_data[q]) : (float) SMALL_LOG_VALUE;

    // Convert to code units.
    my_cloudy.cooling_data[q] -= log10(CoolUnit);
  }
  delete [] temp_data;

  status = H5Dclose(dset_id);
  if (status == h5_error) {
    fprintf(stderr,"Failed to close Cooling dataset.\n");
    return FAIL;
  }

  // Read Heating data.
  if (my_chemistry.UVbackground) {

    temp_data = new double[my_cloudy.data_size];

    sprintf(parameter_name, "/CoolingRates/%s/Heating", group_name);
    dset_id =  H5Dopen(file_id, parameter_name);
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open Heating in %s.\n",my_chemistry.grackle_data_file);
      return FAIL;
    }

    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
    fprintf(stderr,"Reading Cloudy Heating dataset.\n");
    if (status == h5_error) {
      fprintf(stderr,"Failed to read Heating dataset.\n");
      return FAIL;
    }

    my_cloudy.heating_data = new gr_float[my_cloudy.data_size];
    for (q = 0;q < my_cloudy.data_size;q++) {
      my_cloudy.heating_data[q] = temp_data[q] > 0 ? (float) log10(temp_data[q]) : (float) SMALL_LOG_VALUE;

      // Convert to code units.
      my_cloudy.heating_data[q] -= log10(CoolUnit);
    }
    delete [] temp_data;

    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close Heating dataset.\n");
      return FAIL;
    }
  }

  status = H5Fclose (file_id);

  if (my_cloudy.grid_rank > CLOUDY_MAX_DIMENSION) {
    fprintf(stderr,"Error: rank of Cloudy cooling data must be less than or equal to %"ISYM".\n",
	    CLOUDY_MAX_DIMENSION);
    return FAIL;
  }

  return SUCCESS;
}
