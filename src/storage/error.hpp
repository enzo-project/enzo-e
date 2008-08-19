//345678901234567890123456789012345678901234567890123456789012345678901234567890

#ifndef ERROR_HPP
#define ERROR_HPP

/*********************************************************************
 *
 * @file    error.hpp
 * @brief   Error and warning definitions
 * @author  James Bordner 
 * @version 1.0
 *
 * Definitions
 *
 * (*) WARNING(MESSAGE)
 * (*) ERROR(MESSAGE)
 * (*) NOT_IMPLEMENTED(X)
 *
 * $Id$
 *
 */

/// Handle a (non-lethal) warning message

static char warning_message[80];

#define WARNING_MESSAGE(FUNCTION) {		\
    printf ("WARNING File:     %s:%d\n" \
            "WARNING Function: %s()\n" \
            "WARNING Message:  %s\n", \
	    __FILE__,__LINE__,FUNCTION,warning_message); \
  fflush(stdout); \
}

/// Handle a (lethal) error message

static char error_message[80];

#define ERROR_MESSAGE(FUNCTION) { \
    printf ("ERROR File:     %s:%d\n" \
            "ERROR Function: %s()\n" \
            "ERROR Message:  %s\n", \
	    __FILE__,__LINE__,FUNCTION,error_message);			\
  fflush(stdout); \
  exit(1); \
}

/// Placeholder for functions that are not implemented yet

#define NOT_IMPLEMENTED_MESSAGE(X) { \
   printf ("%s:%d WARNING: " X " is not implemented yet\n",__FILE__,__LINE__); \
   fflush(stdout); \
}

#endif
