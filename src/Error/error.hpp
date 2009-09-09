#ifndef ERROR_HPP
#define ERROR_HPP

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

/*********************************************************************
 *
 * @file    error.hpp
 * @brief   Error and warning definitions
 * @author  James Bordner 
 * @version 1.0
 *
 * Definitions
 *
 * (*) WARNING(FUNCTION)
 * (*) ERROR(FUNCTION)
 * (*) INCOMPLETE(FUNCTION)
 *
 * $Id$
 *
 *********************************************************************/

#include "cello.h"
#include "error_exception.hpp"

#define ERROR_MESSAGE_LENGTH 80

//-----------------------------------------------------------------------
/// Handle a (non-lethal) warning message
//-----------------------------------------------------------------------

static char warning_message[ERROR_MESSAGE_LENGTH+1];

#define WARNING_MESSAGE(FUNCTION) {		\
    printf ("\n" \
            "     WARNING File:     %s:%d\n" \
            "     WARNING Function: %s()\n" \
            "     WARNING Message:  %s\n" \
            "\n", \
	    __FILE__,__LINE__,FUNCTION,warning_message); \
  fflush(stdout); \
}

//-----------------------------------------------------------------------
/// Handle a (lethal) error message
//-----------------------------------------------------------------------

static char error_message[ERROR_MESSAGE_LENGTH+1];

#define ERROR_MESSAGE(FUNCTION) { \
    printf ("\n" \
            "     ERROR File:     %s:%d\n" \
            "     ERROR Function: %s()\n" \
            "     ERROR Message:  %s\n" \
            "\n", \
	    __FILE__,__LINE__,FUNCTION,error_message);			\
  fflush(stdout); \
  exit(1); \
}

//-----------------------------------------------------------------------
/// Placeholder for code that is incomplete
//-----------------------------------------------------------------------

static char incomplete_message[ERROR_MESSAGE_LENGTH+1];

#define INCOMPLETE_MESSAGE(FUNCTION) { \
    printf ("\n" \
            "     INCOMPLETE File:     %s:%d\n" \
            "     INCOMPLETE Function: %s()\n" \
            "     INCOMPLETE Message:  %s\n" \
            "\n", \
	    __FILE__,__LINE__,FUNCTION,incomplete_message);			\
   fflush(stdout); \
}

#endif /* ERROR_HPP */

