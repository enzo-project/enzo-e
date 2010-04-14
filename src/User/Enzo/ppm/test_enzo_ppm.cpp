// $Id: test_enzo_method_ppm.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     test_enzo_method_ppm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Apr  1 16:19:18 PDT 2010
/// @brief    Unit tests for the EnzoMethodPpm class

#include <stdio.h>

#include "user.hpp"
#include "test.hpp"

main ()
{
  unit_class ("MethodEnzoPpm");
  MethodEnzoPpm ppm;

  unit_func("initialize");
  ppm.initialize("PPM");
  unit_assert(false);

  unit_func("apply");
  ppm.apply();
  unit_assert(false);

}
