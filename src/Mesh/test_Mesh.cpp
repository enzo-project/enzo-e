// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Mesh class

#include "main.hpp"
#include "test.hpp"
#include "mesh.hpp"

#include "mesh_charm.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class("Mesh");

  unit_func("Mesh");
  Factory * factory = new Factory;
  Mesh * mesh = new Mesh (factory,12,12,12,3,3,3);
  unit_assert(mesh != NULL);

  unit_finalize();

  delete mesh;
  delete factory;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
