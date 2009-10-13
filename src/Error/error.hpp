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

#include <assert.h>

#include "cello.h"
#include "error_exception.hpp"

#define ERROR_MESSAGE_LENGTH 255

//-----------------------------------------------------------------------
/// Handle a (non-lethal) warning message
//-----------------------------------------------------------------------

#define WARNING_MESSAGE(FUNCTION,MESSAGE) { \
    printf ("\n" \
            "     WARNING File:     %s:%d\n" \
            "     WARNING Function: %s()\n" \
            "     WARNING Message:  %s\n" \
            "\n", \
	    __FILE__,__LINE__,FUNCTION,MESSAGE); \
  fflush(stdout); \
}

//-----------------------------------------------------------------------
/// Handle a (lethal) error message
//-----------------------------------------------------------------------

#define ERROR_MESSAGE(FUNCTION,MESSAGE) { \
    printf ("\n" \
            "     ERROR File:     %s:%d\n" \
            "     ERROR Function: %s()\n" \
            "     ERROR Message:  %s\n" \
            "\n", \
	    __FILE__,__LINE__,FUNCTION,MESSAGE);			\
  fflush(stdout); \
  exit(1); \
}

//-----------------------------------------------------------------------
/// Placeholder for code that is incomplete
//-----------------------------------------------------------------------

#define INCOMPLETE_MESSAGE(FUNCTION,MESSAGE) { \
    printf ("\n" \
            "     INCOMPLETE File:     %s:%d\n" \
            "     INCOMPLETE Function: %s()\n" \
            "     INCOMPLETE Message:  %s\n" \
            "\n", \
	    __FILE__,__LINE__,FUNCTION,MESSAGE);			\
   fflush(stdout); \
}

#define TRACE if (trace) {printf ("TRACE %s:%d\n",__FILE__,__LINE__); fflush(stdout); }

#endif /* ERROR_HPP */

