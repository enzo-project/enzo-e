// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ERROR_HPP
#define ERROR_HPP

/// @file     error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-03-18
/// @todo     Add Parallel support
/// @brief    Error- and warning-related defines

#include <assert.h>
#include "cello.h"
#include "error_exception.hpp"

/// @def      ERROR_MESSAGE_LENGTH
/// @brief    Maximum length of error and warning messages
#define ERROR_MESSAGE_LENGTH 255

/// @def      WARNING_MESSAGE
/// @brief    Handle a (non-lethal) warning message
#define WARNING_MESSAGE(FUNCTION,MESSAGE) {		\
    printf ("\n"					\
            "     WARNING File:     %s:%d\n"		\
            "     WARNING Function: %s()\n"		\
            "     WARNING Message:  %s\n"		\
            "\n",					\
	    __FILE__,__LINE__,FUNCTION,MESSAGE);	\
    fflush(stdout);					\
  }

/// @def      ERROR_MESSAGE
/// @brief    Handle a (lethal) error message
#define ERROR_MESSAGE(FUNCTION,MESSAGE) {		\
    printf ("\n"					\
            "     ERROR File:     %s:%d\n"		\
            "     ERROR Function: %s()\n"		\
            "     ERROR Message:  %s\n"			\
            "\n",					\
	    __FILE__,__LINE__,FUNCTION,MESSAGE);	\
    fflush(stdout);					\
    exit(1);						\
  }

/// @def      INCOMPLETE_MESSAGE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE_MESSAGE(FUNCTION,MESSAGE) {		\
    printf ("\n"					\
            "     INCOMPLETE File:     %s:%d\n"		\
            "     INCOMPLETE Function: %s()\n"		\
            "     INCOMPLETE Message:  %s\n"		\
            "\n",					\
	    __FILE__,__LINE__,FUNCTION,MESSAGE);	\
    fflush(stdout);					\
  }

/// @def      TRACE
/// @brief    Trace file name and location to stdout
#define TRACE {								\
    if (trace) {							\
      printf ("TRACE %s:%d\n",__FILE__,__LINE__); fflush(stdout);	\
    }									\
  }

#endif /* ERROR_HPP */

