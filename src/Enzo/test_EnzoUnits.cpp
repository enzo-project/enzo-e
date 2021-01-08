// See LICENSE_CELLO file for license and copyright information

/// @file     test_EnzoUnits.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-12
/// @brief    Test program for the EnzoUnits class

#include "test.hpp"
#include "main.hpp"
#include "enzo.hpp"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class ("EnzoUnits");

  EnzoUnits * units = new EnzoUnits;

  unit_assert (units != NULL);

  unit_class ("EnzoPhysicsCosmology");

  EnzoPhysicsCosmology * cosmology = new EnzoPhysicsCosmology;

  unit_assert (cosmology != NULL);

  //  units->set_cosmology(cosmology);
 
  units->set_using_mass(1.0,4.0,8.0);

  unit_assert (units->length() == 1.0);
  unit_assert (units->mass() == 4.0);
  unit_assert (units->time() == 8.0);

  enzo_float omega_baryon_now  =  0.04;
  enzo_float omega_cdm_now     =  0.26;
  
  cosmology->set_comoving_box_size ( 64.0);
  cosmology->set_omega_matter_now    ( 0.3);
  cosmology->set_omega_lambda_now    ( 0.7);
  cosmology->set_hubble_constant_now ( 0.5);
  cosmology->set_max_expansion_rate  ( 0.015);
  cosmology->set_initial_redshift    (30.0);
  cosmology->set_final_redshift      ( 0.0);

  enzo_float t0 = cosmology->time_from_redshift(cosmology->initial_redshift());
  unit_func ("time_from_redshift()");
  unit_assert (cello::err_rel(t0,(enzo_float)0.81650250388244) <= 1e-15);
  cosmology->set_current_time(0.0);

  cosmology->print();

  CkPrintf ("length = %g\n",units->length());
  CkPrintf ("mass   = %g\n",units->mass());
  CkPrintf ("time   = %g\n",units->time());
  
  CkPrintf ("volume = %g\n",units->volume());
  
  unit_finalize();
  
  exit_();
}


PARALLEL_MAIN_END
#include "enzo.def.h"

