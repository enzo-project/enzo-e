// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Mesh class

#include "test.hpp"

#include "mesh.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Mesh.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();
  unit_class ("Mesh");

  DataDescr data_descr;

  unit_func("Mesh");
  Mesh * mesh = new Mesh (&data_descr);
  unit_assert(mesh != NULL);

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Mesh.def.h)
