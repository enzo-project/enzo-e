#ifndef TEST_UNIT_HPP
#define TEST_UNIT_HPP

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
