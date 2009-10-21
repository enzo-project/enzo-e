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
 
#include "cello.h"

#include "error.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{
  bool passed = false;

  unit_class ("Control");
  unit_open();

  //  Control control;

  unit_assert(passed);

  unit_close();
}
