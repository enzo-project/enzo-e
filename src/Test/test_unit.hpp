/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * SYNOPSIS:
 *
 *    
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

#ifndef TEST_UNIT_HPP
#define TEST_UNIT_HPP

/**
 *********************************************************************
 *
 * @file    unit.hpp
 * @brief   Unit testing functions
 * @author  James Bordner 
 * @date    Sat Feb 23 15:22:59 PST 2008
 *
 * $Id$
 *
 *********************************************************************
 */

//======================================================================
// INCLUDES
//======================================================================

#include <string.h>

#include "error.hpp"

//======================================================================
// DEFINES
//======================================================================


//----------------------------------------------------------------------

#define UNIT_MAX_NAME_LEN 40

//======================================================================
// VARIABLES
//======================================================================

namespace unit {

  char class_name[UNIT_MAX_NAME_LEN] = {0};
  char func_name[UNIT_MAX_NAME_LEN] = {0};

  FILE *fp = 0;

  int test_num = 1;

  const char * pass_string = "\033[01;32mPass\033[00m";
  const char * fail_string = "\033[01;31mFAIL\033[00m";
  //const char * pass_string = " [ Pass ]";
  //const char * fail_string = " [ FAIL ]";

}

//======================================================================
// FUNCTIONS
//======================================================================

void unit_class (const char * c)
{
  strncpy (unit::class_name,c,UNIT_MAX_NAME_LEN);
}

//----------------------------------------------------------------------

void unit_func (const char * f)
{
  strncpy (unit::func_name,f,UNIT_MAX_NAME_LEN);
}

//----------------------------------------------------------------------


#define unit_assert(RESULT) unit_assert_(RESULT, __FILE__,__LINE__);

void unit_assert_ (bool result, const char * file, int line)
{
  printf ("%s %s:%3d  %s::%s() %d\n",
	  (result)? unit::pass_string : unit::fail_string,
	  file,line,unit::class_name,unit::func_name,unit::test_num);
  if (unit::fp) {
    fprintf (unit::fp,"%s %s:%3d  %s::%s() %d\n",
	     (result)? unit::pass_string : unit::fail_string,
	     file,line,unit::class_name,unit::func_name,unit::test_num);
  } else {
    WARNING_MESSAGE("Test::unit_assert","unit_open() not called before unit_assert()"); 
  }
  unit::test_num++;
}

//----------------------------------------------------------------------

void unit_open ()
{
  char filename [UNIT_MAX_NAME_LEN+5];
  if (strlen(unit::class_name)==0) {
    WARNING_MESSAGE("Test::unit_open","unit_class() not called before unit_open()"); 
  }
  sprintf (filename,"%s.unit",unit::class_name);
  unit::fp = fopen (filename,"w");
}

//----------------------------------------------------------------------

void unit_close ()
{
  if (unit::fp) {
    fclose (unit::fp);
  } else {
    WARNING_MESSAGE("Test::unit_close","unit_open() not called before unit_close()"); 
  }
}

#endif
