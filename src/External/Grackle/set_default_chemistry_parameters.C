/***********************************************************************
/
/ Set default parameter values
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "chemistry_data.h"
#include "code_units.h"

chemistry_data set_default_chemistry_parameters()
{
  
  chemistry_data my_chemistry;

  my_chemistry.Gamma                          = 5./3.;
  my_chemistry.use_grackle                    = FALSE;  // off
  my_chemistry.with_radiative_cooling         = TRUE;
  my_chemistry.primordial_chemistry           = FALSE;  // off
  my_chemistry.metal_cooling                  = FALSE;
  my_chemistry.h2_on_dust                     = FALSE;

  my_chemistry.cmb_temperature_floor          = TRUE;   // use CMB floor.
  my_chemistry.grackle_data_file              = "";

  my_chemistry.three_body_rate                = 0;   // ABN02
  my_chemistry.cie_cooling                    = 1;
  my_chemistry.h2_optical_depth_approximation = 1;

  my_chemistry.photoelectric_heating          = 0;
  my_chemistry.photoelectric_heating_rate     = 8.5e-26;  // ergs cm-3 s-1

  my_chemistry.UVbackground                   = 0;

  my_chemistry.UVbackground_table.Nz     = 0;
  my_chemistry.UVbackground_table.z      = NULL;
  my_chemistry.UVbackground_table.k24    = NULL;
  my_chemistry.UVbackground_table.k25    = NULL;
  my_chemistry.UVbackground_table.k26    = NULL;
  my_chemistry.UVbackground_table.k27    = NULL;
  my_chemistry.UVbackground_table.k28    = NULL;
  my_chemistry.UVbackground_table.k29    = NULL;
  my_chemistry.UVbackground_table.k30    = NULL;
  my_chemistry.UVbackground_table.k31    = NULL;
  my_chemistry.UVbackground_table.piHI   = NULL;
  my_chemistry.UVbackground_table.piHeII = NULL;
  my_chemistry.UVbackground_table.piHeI  = NULL;

  my_chemistry.UVbackground_redshift_on      = 7.0;
  my_chemistry.UVbackground_redshift_off     = 0.0;
  my_chemistry.UVbackground_redshift_fullon  = 6.0;
  my_chemistry.UVbackground_redshift_drop    = 0.0;

  my_chemistry.Compton_xray_heating   = 0;

  my_chemistry.LWbackground_intensity = 0.0;   // [in units of 10^21 erg/s/cm^2/Hz/sr]
  my_chemistry.LWbackground_sawtooth_suppression = 0;

  my_chemistry.HydrogenFractionByMass       = 0.76;
  /* The DToHRatio is by mass in the code, so multiply by 2. */
  my_chemistry.DeuteriumToHydrogenRatio     = 2.0*3.4e-5; // Burles & Tytler 1998
  my_chemistry.SolarMetalFractionByMass     = 0.02041;
  my_chemistry.NumberOfTemperatureBins      = 600;
  my_chemistry.ih2co                        = 1;
  my_chemistry.ipiht                        = 1;
  my_chemistry.TemperatureStart             = 1.0;
  my_chemistry.TemperatureEnd               = 1.0e9;
  my_chemistry.comp_xray                    = 0;
  my_chemistry.temp_xray                    = 0;
  my_chemistry.CaseBRecombination           = 0; // default to case A rates
  my_chemistry.NumberOfDustTemperatureBins  = 250;
  my_chemistry.DustTemperatureStart         = 1.0;
  my_chemistry.DustTemperatureEnd           = 1500.0;

  my_chemistry.cloudy_metal.grid_rank        = 0;
  my_chemistry.cloudy_electron_fraction_factor = 9.153959e-3; // Cloudy 07.02 abundances

  return my_chemistry;
}
