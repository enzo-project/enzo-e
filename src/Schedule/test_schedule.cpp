// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      test_schedule.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:47:35 PST 2008
/// @brief     Program implementing unit tests for the Schedule class
 
#include "cello.h"

#include "error.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{
  bool passed = false;

  unit_class ("Schedule");
  unit_func ("null");
  unit_open();

  //  Schedule schedule;

  //FAILS
  unit_assert(passed);

  unit_close();
}
