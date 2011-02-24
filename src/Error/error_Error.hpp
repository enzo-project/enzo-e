// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ERROR_ERROR_HPP
#define ERROR_ERROR_HPP

/// @file     error_Error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Add Parallel support
/// @bug      exit() is called instead of MPI_Abort(), etc.
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the Error class
/// @def      ERROR_LENGTH
/// @brief    Maximum length of error and warning messages
#define ERROR_LENGTH 255

/// @def      WARNING
/// @brief    Handle a (non-lethal) warning message
#define WARNING(FUNCTION,MESSAGE)					\
  {									\
    message_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
  }

/// @def      UNTESTED
/// @brief    Handle a (non-lethal) untested message
#define UNTESTED(FUNCTION)					\
  {								\
    message_(stdout,"UNTESTED",__FILE__,__LINE__,FUNCTION,"");	\
  }

/// @def      ERROR
/// @brief    Handle a (lethal) error message
#define ERROR(FUNCTION,MESSAGE)						\
  {									\
    message_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
    exit(1);								\
  }

/// @def      INCOMPLETE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE(FUNCTION,MESSAGE)					\
  {									\
    message_(stdout,"INCOMPLETE",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
  }

/// @def      TRACE
/// @brief    Trace file name and location to stdout
#define TRACE(MESSAGE)						\
  {								\
    PARALLEL_PRINTF ("TRACE %s:%d %s\n",__FILE__,__LINE__,MESSAGE);	\
    fflush(stdout);						\
  }

/// @def      ASSERT
/// @brief    Equivalent to assert()
#define ASSERT(FUNCTION,MESSAGE,ASSERTION)				\
  {									\
    if (!(ASSERTION)) {							\
      message_(stderr,"ASSERT",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
      PARALLEL_EXIT;							\
    }									\
  }


#define message_(FP,TYPE,FILE,LINE,FUNCTION,MESSAGE)		\
  {								\
    fprintf (FP,"     ======================================================================\n");						\
    fprintf (FP,"     %s File:     %s:%d\n",TYPE,FILE,LINE);	\
    if (strcmp(FUNCTION,"") != 0)				\
      fprintf (FP,"     %s Function: %s()\n", TYPE,FUNCTION);	\
    if (strcmp(MESSAGE,"") != 0)				\
      fprintf (FP,"     %s Message:  %s\n",   TYPE,MESSAGE);	\
    fprintf (FP,"     ======================================================================\n");						\
  }

#endif /* ERROR_ERROR_HPP */

