// See LICENSE_CELLO file for license and copyright information

/// @file     test_It.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the It class

#include "main.hpp"
#include "test.hpp"

#include "component.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  GroupProcess * group_process = GroupProcess::create();

  unit_init(0,1);

  // Create Hierarchy over which to iterate

  Factory * factory = new Factory;
  Hierarchy * hierarchy = new Hierarchy (factory);
  unit_assert(hierarchy != NULL);

  FieldDescr field_descr;
  GroupProcess * group_process = GroupProcess::create();
  hierarchy->create_root_patch(group_process,&field_descr,12,12,12,3,3,3);
  unit_assert(hierarchy->patch(0)!=NULL);


  FieldDescr * field_descr = new FieldDescr;

  // Set Patch size (12,12,12)

  int patch_size[] = {12,12,12};

  int patch_blocking[] = {3,3,3};

  // Set domain extents

  double domain_lower[] = {0.0, 0.0, 0.0};
  double domain_upper[] = {1.0, 1.0, 1.0};

  Patch * patch = factory->create_patch 
    (group_process,
     patch_size[0],     patch_size[1],     patch_size[2],
     patch_blocking[0], patch_blocking[1], patch_blocking[2],
     domain_lower[0],   domain_lower[1],   domain_lower[2],
     domain_upper[0],   domain_upper[1],   domain_upper[2]);

  unit_class("It");

  ItPatch * it_patch = new ItPatch;

  unit_func ("function");

  unit_assert (it != NULL)

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

