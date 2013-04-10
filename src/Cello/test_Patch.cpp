// See LICENSE_CELLO file for license and copyright information

/// @file     test_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Patch class

#ifdef REMOVE_PATCH
#else /* REMOVE_PATCH */


#include <assert.h>
#include "main.hpp" 
#include "test.hpp"

#include "mesh.hpp"

#ifdef CONFIG_USE_CHARM
CProxy_Patch * patch_global;
#else
Patch * patch_global;
#endif

#ifdef CONFIG_USE_CHARM
extern CProxy_SimulationCharm proxy_simulation;
#endif

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Patch");

  //======================================================================
  // np = 1
  // ip = 0 [default]
  //======================================================================

  unit_func("Patch");

  const FieldDescr * field_descr = new FieldDescr;

  Factory * factory = new Factory;

  // Set Patch size, offset, and blocking

  const int patch_size[]     = {12,12,12};
  const int patch_offset[]   = {5, 2, 9};
  const int patch_blocking[] = {3,3,3};

  // Set domain extents

  const double domain_lower[] = {0.0, 0.0, 0.0};
  const double domain_upper[] = {1.0, 1.0, 1.0};

  patch_global = factory->create_patch 
    (
     field_descr,
     patch_size[0],     patch_size[1],     patch_size[2],
     patch_offset[0],   patch_offset[1],   patch_offset[2],
     patch_blocking[0], patch_blocking[1], patch_blocking[2],
     domain_lower[0],   domain_lower[1],   domain_lower[2],
     domain_upper[0],   domain_upper[1],   domain_upper[2]);

  patch_global->p_test();

  delete factory;
  delete field_descr;
}

PARALLEL_MAIN_END

#endif /* REMOVE_PATCH */
