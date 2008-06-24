#ifndef UNIT_HPP
#define UNIT_HPP

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

//======================================================================
// DEFINES
//======================================================================


#define unit_assert(RESULT)		\
  printf ("%s %s:%3d  %s::%s() %d\n",	\
	  (RESULT)? unit::pass_string : unit::fail_string,			\
	  __FILE__,__LINE__,unit::class_name,unit::func_name,unit::test_num); \
  fprintf (unit::fp,"%s %s:%3d  %s::%s() %d\n",				\
	  (RESULT)? unit::pass_string : unit::fail_string,			\
	  __FILE__,__LINE__,unit::class_name,unit::func_name,unit::test_num); \
  unit::test_num++;

//----------------------------------------------------------------------

#define unit_class_size(CLASS) \
  printf ("sizeof(%s) = %d\n",unit::class_name,sizeof(CLASS)); \
  fprintf (unit::fp,"sizeof(%s) = %d\n",unit::class_name,sizeof(CLASS));


#define UNIT_MAX_NAME_LEN 40

//======================================================================
// VARIABLES
//======================================================================

namespace unit {

  char class_name[UNIT_MAX_NAME_LEN];
  char func_name[UNIT_MAX_NAME_LEN];

  FILE *fp;

  int test_num = 1;

  const char * pass_string = "\033[01;32mPass\033[00m";
  const char * fail_string = "\033[01;38mFail\033[00m";

}

//======================================================================
// FUNCTIONS
//======================================================================

void unit_class (char * c)
{
  strncpy (unit::class_name,c,UNIT_MAX_NAME_LEN);
}

//----------------------------------------------------------------------

void unit_func (char * f)
{
  strncpy (unit::func_name,f,UNIT_MAX_NAME_LEN);
}

//----------------------------------------------------------------------

void unit_open ()
{
  char filename [UNIT_MAX_NAME_LEN+5];
  sprintf (filename,"%s.unit",unit::class_name);
  unit::fp = fopen (filename,"w");
}

//----------------------------------------------------------------------

void unit_close ()
{
  fclose (unit::fp);
}

#endif
