/** 
 *********************************************************************
 *
 * @file      test_perf.cpp
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

#include "test_unit_.hpp"
#include "data_scalar_.hpp"
#include "data_array_.hpp"
#include "timer.hpp"

main()
{
  Timer timer;
  // Timer tests
  unit_class ("Timer");
  unit_open();
  unit_class_size(Timer);

  printf ("Initial timer value = %24.16f\n",timer.value());

  timer.start();
  system("sleep 1");
  timer.stop();

  printf ("Initial timer value = %24.16f\n",timer.value());

  unit_assert((timer.value() - 1.0) < 0.01);

  timer.start();
  system("sleep 1");
  timer.stop();

  printf ("Initial timer value = %24.16f\n",timer.value());

  unit_assert((timer.value() - 2.0) < 0.01);

  unit_close();
}
