/***********************************************************************
/
/ Code units data structure
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __CODE_UNITS_H__
#define __CODE_UNITS_H__
typedef struct
{
  gr_int comoving_coordinates;
  gr_float density_units;
  gr_float length_units;
  gr_float time_units;
  gr_float velocity_units;
  gr_float a_units;
} code_units;
#endif
