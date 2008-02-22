//345678901234567890123456789012345678901234567890123456789012345678901234567890

#ifndef UNIT_HPP
#define UNIT_HPP

// $Id$
/**
 * @file    unit.hpp
 * @brief   Unit testing functions
 * @author  James Bordner 
 * @version 1.0
 *
 * Defines
 *
 * ( ) UNIT_ASSERT(RESULT)
 * ( ) UNIT_CLASS (char * class)
 * ( ) UNIT_FUNC (char * func)
 *
 */
// $Log$

#include <string.h>

#define UNIT_ASSERT(RESULT)		\
  printf ("%s %s:%3d  %s::%s() %d\n",	\
	  (RESULT)?" PASS ":" FAIL ",__FILE__,__LINE__,class_name,func_name,test_num++);


#define UNIT_MAX_NAME_LEN 40

char class_name[UNIT_MAX_NAME_LEN];
char func_name[UNIT_MAX_NAME_LEN];

int test_num = 1;

void UNIT_CLASS (char * c)
{
  strncpy (class_name,c,UNIT_MAX_NAME_LEN);
};

void UNIT_FUNC (char * f)
{
  strncpy (func_name,f,UNIT_MAX_NAME_LEN);
};

#endif
