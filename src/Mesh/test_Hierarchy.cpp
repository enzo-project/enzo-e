// See LICENSE_CELLO file for license and copyright information

/// @file     test_Hierarchy.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Hierarchy class

#include "main.hpp"
#include "test.hpp"
#include "mesh.hpp"

#include "mesh_charm.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;
  unit_init();

  unit_class("Hierarchy");

  unit_func("Hierarchy");
  Factory * factory = new Factory;
  Hierarchy * hierarchy = new Hierarchy (factory);
  unit_assert(hierarchy != NULL);

  FieldDescr field_descr;
  GroupProcess * group_process = GroupProcess::create();
  hierarchy->create_root_patch(group_process,&field_descr,12,12,12,3,3,3);
  unit_assert(hierarchy->patch(0)!=NULL);

    
  unit_finalize();

  delete hierarchy;
  delete factory;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
