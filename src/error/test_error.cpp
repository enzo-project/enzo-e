/** 
 *********************************************************************
 *
 * @file      test_error.cpp
 * @brief     Program implementing unit tests for error classes
 * @author    James Bordner
 * @date      Wed Aug 20 11:24:14 PDT 2008
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "error.hpp"
#include "test_unit_.hpp"

main()
{

  unit_class ("Error");
  unit_open();

  //----------------------------------------------------------------------
  printf ("Warning message:\n");

  sprintf (warning_message,"Warning message test");
  WARNING_MESSAGE("main");

  unit_assert (true);

  //----------------------------------------------------------------------
  printf ("Incomplete message:\n");

  sprintf (incomplete_message,"Incomplete message test");
  INCOMPLETE_MESSAGE("main");

  unit_assert (true);

  //----------------------------------------------------------------------
  printf ("Error message:\n");

  sprintf (error_message,"Error message test");
  ERROR_MESSAGE("main");
  
  // Errors should abort, so all following linse should not be executed
  unit_assert (false); 
  //----------------------------------------------------------------------

  unit_close();
}
