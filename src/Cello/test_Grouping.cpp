// See LICENSE_CELLO file for license and copyright information

/// @file     test_Grouping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 15:35:51 PDT 2014
/// @brief    Test program for the Grouping class

#include "main.hpp"
#include "test.hpp"
#include <stdexcept>
#include <string>
#include "field_Grouping.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Grouping");

  Grouping * grouping = new Grouping;

  unit_func("Grouping::size()");
  unit_assert(grouping->size("group_1")==0);

  unit_func("Grouping::add()");
  grouping->add("density",   "group_1");
  unit_func("Grouping::size()");

  unit_func("Grouping::add()");
  grouping->add("velocity_x","group_2");
  unit_func("Grouping::add()");
  grouping->add("velocity_y","group_2");
  grouping->add("velocity_y","group_2");
  grouping->add("velocity_y","group_2");
  unit_func("Grouping::add()");
  grouping->add("velocity_z","group_2");

  unit_func("Grouping::add()");
  grouping->add("density",   "group_3");
  unit_func("Grouping::add()");
  grouping->add("velocity_x","group_3");
  unit_func("Grouping::add()");
  grouping->add("velocity_y","group_3");
  unit_func("Grouping::add()");
  grouping->add("velocity_z","group_3");

  unit_func("Grouping::size()");
  unit_assert(grouping->size("group_1")==1);
  unit_assert(grouping->size("group_2")==3);
  unit_assert(grouping->size("group_3")==4);

  unit_func("Grouping::is_in()");

  // non-field tests
  unit_assert(grouping->is_in("not_a_field", "group_1") == false);
  unit_assert(grouping->is_in("not_a_field", "group_2") == false);

  // non-group tests
  unit_assert(grouping->is_in("density",    "not_a_group") == false);
  unit_assert(grouping->is_in("velocity_x", "not_a_group") == false);
  unit_assert(grouping->is_in("velocity_y", "not_a_group") == false);
  unit_assert(grouping->is_in("velocity_z", "not_a_group") == false);

  // non-field and non-group
  unit_assert(grouping->is_in("not_a_field", "not_a_group") == false);
  unit_assert(grouping->is_in("not_a_field", "not_a_group") == false);


  // group 1 tests
  unit_assert(grouping->is_in("density",   "group_1") == true);
  unit_assert(grouping->is_in("velocity_x","group_1") == false);
  unit_assert(grouping->is_in("velocity_y","group_1") == false);
  unit_assert(grouping->is_in("velocity_z","group_1") == false);

  // group 2 tests
  unit_assert(grouping->is_in("density",   "group_2") == false);
  unit_assert(grouping->is_in("velocity_x","group_2") == true);
  unit_assert(grouping->is_in("velocity_y","group_2") == true);
  unit_assert(grouping->is_in("velocity_z","group_2") == true);

  // group 3 tests
  unit_assert(grouping->is_in("density",   "group_3") == true);
  unit_assert(grouping->is_in("velocity_x","group_3") == true);
  unit_assert(grouping->is_in("velocity_y","group_3") == true);
  unit_assert(grouping->is_in("velocity_z","group_3") == true);

  unit_func("field_in_group()");

  std::string g1_0 = grouping->item("group_1",0);
  std::string g1_1 = grouping->item("group_1",1);
  unit_assert(g1_0 == "density");
  unit_assert(g1_1 == "");
  

  //--------------------------------------------------

  delete grouping;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
