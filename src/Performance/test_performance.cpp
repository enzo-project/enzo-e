//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      test_performance.cpp
 * @brief     Program implementing unit tests for performance classes
 * @author    James Bordner
 * @date      Wed Apr 23 12:25:18 PDT 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "test.hpp"
#include "array.hpp"
#include "performance.hpp"

int main(int argc, char ** argv)
{
  Timer timer;
  const double time_tolerance = 0.05;

  // Timer tests

  unit_class ("Timer");
  unit_open();

  printf ("Initial timer value = %24.16f\n",timer.value());

  int return_value;
  timer.start();

  return_value = system("sleep 1");
  timer.stop();

  printf ("Initial timer value = %24.16f\n",timer.value());

  unit_assert((timer.value() - 1.0) < time_tolerance);

  timer.start();
  return_value = system("sleep 1");
  timer.stop();

  printf ("Initial timer value = %24.16f\n",timer.value());

  unit_assert((timer.value() - 2.0) < time_tolerance);

  unit_close();
}
