// $Id: test_enzo_method_ppm.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_enzo_method_ppm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Apr  1 16:19:18 PDT 2010
/// @brief    Unit tests for the EnzoMethodPpm class

#include <stdio.h>

#include "method.hpp"
#include "test.hpp"

main ()
{
  unit_class ("EnzoMethodPpm");
  unit_open();
  unit_func("initialize");
  unit_assert(false);
  unit_close();
}
