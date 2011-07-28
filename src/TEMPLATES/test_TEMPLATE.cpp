// See LICENSE_CELLO file for license and copyright information

/// @file     test_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Classname class

#include "test.hpp"

#include "component.hpp"

#include PARALLEL_CHARM_INCLUDE(test.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class("Classname");

  Classname * classname = new Classname;

  unit_func ("function");

  unit_assert (classname != NULL)

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test.def.h)
