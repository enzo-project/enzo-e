// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef TEST_UNIT_HPP
#define TEST_UNIT_HPP

/// @file     test_unit.hpp
/// @brief    Define the unit namespace and unit testing functions
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Convert namespace to class
/// @todo     Address need to edit code to use pass/fail string colors or not
/// @date     Sat Feb 23 15:22:59 PST 2008

#include <string.h>
#include "error.hpp"

/// @def      UNIT_MAX_NAME_LEN
/// @brief    Maximum length of a class or function name
#define UNIT_MAX_NAME_LEN 40

namespace unit {

  /// @namespace unit
  /// @ingroup   Test
  /// @brief     Current class name, function name, and test results
       
  char class_name[UNIT_MAX_NAME_LEN] = {0};
  char func_name[UNIT_MAX_NAME_LEN] = {0};

  FILE *fp = 0;

  int test_num = 1;

  const char * pass_string = "\033[01;32mPass\033[00m";
  const char * fail_string = "\033[01;31mFAIL\033[00m";
  //const char * pass_string = " [ Pass ]";
  //const char * fail_string = " [ FAIL ]";

}

/// @brief Set the current unit testing class name
void unit_class (const char * c)
{
  strncpy (unit::class_name,c,UNIT_MAX_NAME_LEN);
}

/// @brief Set the current unit testing function name
void unit_func (const char * f)
{
  strncpy (unit::func_name,f,UNIT_MAX_NAME_LEN);
}

/// @brief Assert result of test; macro used for FILE and LINE expansion
#define unit_assert(RESULT) unit_assert_(RESULT, __FILE__,__LINE__);

/// @brief Assert result of test macro; called by unit_assert macro
void unit_assert_ (bool result, const char * file, int line)
{
  printf ("%s %s:%d  %s::%s() %d\n",
	  (result)? unit::pass_string : unit::fail_string,
	  file,line,unit::class_name,unit::func_name,unit::test_num);
  if (unit::fp) {
    fprintf (unit::fp,"%s %s:%d  %s::%s() %d\n",
	     (result)? unit::pass_string : unit::fail_string,
	     file,line,unit::class_name,unit::func_name,unit::test_num);
  } else {
    WARNING_MESSAGE("Test::unit_assert","unit_open() not called before unit_assert()"); 
  }
  unit::test_num++;
}

/// @brief Open the file for unit test results.  Call to unit_class() required
void unit_open ()
{
  char filename [UNIT_MAX_NAME_LEN+5];
  if (strlen(unit::class_name)==0) {
    WARNING_MESSAGE("Test::unit_open","unit_class() not called before unit_open()"); 
  }
  sprintf (filename,"%s.unit",unit::class_name);
  unit::fp = fopen (filename,"w");
}

/// @brief Close the file for unit test results.  Call to unit_open() required
void unit_close ()
{
  if (unit::fp) {
    fclose (unit::fp);
  } else {
    WARNING_MESSAGE("Test::unit_close","unit_open() not called before unit_close()"); 
  }
}

#endif
