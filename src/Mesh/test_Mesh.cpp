// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Mesh class
 
#include <stdio.h>
#include <string>

#include "mesh.hpp"

#include "test.hpp"
#include "parallel.def"

#include PARALLEL_CHARM_INCLUDE(test_Mesh.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();
  unit_class ("Mesh");
  unit_func("Mesh");

  DataDescr data_descr;

  Mesh mesh(&data_descr);

  unit_assert(false);
  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Mesh.def.h)
