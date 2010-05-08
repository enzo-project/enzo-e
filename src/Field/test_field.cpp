// $Id: test_array.cpp 1412 2010-05-05 00:10:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_field.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Test program for the FieldDescr class

#include "cello.h"

#include "test.hpp"
#include "field.hpp"

int main()
{

  unit_class ("Field");
  unit_func ("Field");
  FieldDescr * field = new FieldDescr;
  unit_assert(field != 0);

}
