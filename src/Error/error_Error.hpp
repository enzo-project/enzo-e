// See LICENSE_CELLO file for license and copyright information

#ifndef ERROR_ERROR_HPP
#define ERROR_ERROR_HPP

//----------------------------------------------------------------------
/// @file     error_Error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      exit() is called instead of MPI_Abort(), etc.
/// @date     2011-04-07
/// @brief    Declaration of the Error class
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// @def      ERROR_LENGTH
/// @brief    Maximum length of error and warning messages
#define ERROR_LENGTH 255

//----------------------------------------------------------------------
/// @def      WARNING
/// @brief    Handle a (non-lethal) warning message
#define WARNING(FUNCTION,MESSAGE)					\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
  }

//----------------------------------------------------------------------
/// @def      UNTESTED
/// @brief    Handle a (non-lethal) untested message
#define UNTESTED(FUNCTION)					\
  {								\
    message2_(stdout,"UNTESTED",__FILE__,__LINE__,FUNCTION,"");	\
  }

//----------------------------------------------------------------------
/// @def      ERROR
/// @brief    Handle a (lethal) error message
#define ERROR(FUNCTION,MESSAGE)						\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
    exit(1);								\
  }
#define ERROR2(FUNCTION,MESSAGE,ARG1)					\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,ARG1);	\
    exit(1);								\
  }

//----------------------------------------------------------------------
/// @def      INCOMPLETE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE(FUNCTION)					\
  {									\
    message2_(stdout,"INCOMPLETE",__FILE__,__LINE__,FUNCTION,"");	\
  }

//----------------------------------------------------------------------
/// @def      TRACE
/// @brief    Trace file name and location to stdout
#ifdef CELLO_TRACE
#   define TRACE(MESSAGE)						\
  {									\
    PARALLEL_PRINTF ("TRACE %s:%d %s\n",__FILE__,__LINE__,MESSAGE);		\
    fflush(stdout);						\
  }
#else
#   define TRACE(MESSAGE) /* This space intentionally left blank */
#endif

//----------------------------------------------------------------------
/// @def      ASSERT
/// @brief    Equivalent to assert()
#define ASSERT(FUNCTION,MESSAGE,ASSERTION)				\
  {									\
    if (!(ASSERTION)) {							\
      message2_(stderr,"ASSERT",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
      exit(1); /* VIOLATES PARALLEL_EXIT() */                           \
    }									\
  }


//----------------------------------------------------------------------
/// @def message_
/// @brief write the given error, warning, etc. message
#define message_(FP,TYPE,FILE,LINE,FUNCTION,MESSAGE)		\
  {								\
    fprintf (FP,"\n");			\
    fprintf (FP,"     %10s  %s:%d\n",TYPE,FILE,LINE);	\
    if (strcmp(FUNCTION,"") != 0)				\
      fprintf (FP,"     %10s  %s()\n", TYPE,FUNCTION);	\
    if (strcmp(MESSAGE,"") != 0)				\
      fprintf (FP,"     %10s  %s\n",   TYPE,MESSAGE);	\
    fprintf (FP,"\n");						\
  }

extern void message2_
(FILE * fp,
 const char * type, 
 const char * file, 
 int line, 
 const char * function, 
 const char * message,
 ...);

#endif /* ERROR_ERROR_HPP */

