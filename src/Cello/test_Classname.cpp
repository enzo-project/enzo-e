// See LICENSE_CELLO file for license and copyright information

/// @file     test_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Classname class

#include "main.hpp"
#include "test.hpp"

#include "component.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Classname");

  Classname * classname = new Classname;

  unit_assert (classname != NULL);

  //--------------------------------------------------

  unit_func ("function()");

  unit_assert (false);

  //--------------------------------------------------

  delete classname;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

