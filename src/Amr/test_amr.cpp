// $Id: test_amr.cpp 1302 2010-03-17 00:16:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_amr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Amr class
 
#include <stdio.h>
#include <string>

#include "cello.h"
#include "test.hpp"
#include "amr.hpp"

int main(int argc, char ** argv)
{
  unit_class ("Amr");
  unit_func("Amr");
  Amr amr;
  unit_assert(false);
}
