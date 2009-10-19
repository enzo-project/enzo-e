//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      test_control.cpp
 * @brief     Program implementing unit tests for the Control class
 * @author    James Bordner
 * @date      Thu Feb 21 16:47:35 PST 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <string>

#include "error.hpp"
#include "test.hpp"
// #include "control.hpp"

int main(int argc, char ** argv)
{
  bool passed = false;

  unit_class ("Control");
  unit_open();

  //  Control control;

  unit_assert(passed);

  unit_close();
}
