//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      test_parameters.cpp
 * @brief     Program implementing unit tests for the Parameters class
 * @author    James Bordner
 * @date      Thu Feb 21 16:04:03 PST 2008
 * 
 * $Id: test_parameters.cpp 715 2009-07-08 23:48:09Z bordner $
 * 
 *********************************************************************
 */
 
#include <stdio.h>

#include "test.hpp"
#include "error.hpp"
#include "parameters.hpp"

main()
{

  unit_class ("Parameters");

  unit_open();

  //----------------------------------------------------------------------
  // test parameter
  //----------------------------------------------------------------------

  Parameters * parameters = new Parameters;

  unit_func("set_group()");

  printf ("parameters->get_group() = %s\n",parameters->get_group().c_str());
  parameters->set_group("Group 1");

  unit_assert(parameters->get_group() == "Group 1");
  unit_assert(parameters->get_subgroup() == "");

  unit_func("read_bison()");
  FILE *fp = fopen ("implosion-1.0.in","r");
  parameters->read ( fp );

  unit_close();

  delete parameters;

}
