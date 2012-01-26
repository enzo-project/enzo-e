// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     error_Error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-07
/// @brief    Declaration of the Error class
//----------------------------------------------------------------------

#ifndef ERROR_ERROR_HPP
#define ERROR_ERROR_HPP

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

#define WARNING1(FUNCTION,MESSAGE,ARG1)					\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1);							\
  }

#define WARNING2(FUNCTION,MESSAGE,ARG1,ARG2)				\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2);						\
  }
#define WARNING3(FUNCTION,MESSAGE,ARG1,ARG2,ARG3)			\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3);						\
  }
#define WARNING4(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4)			\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4);					\
  }
#define WARNING5(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5)		\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4,ARG5);				\
  }
#define WARNING6(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6)	\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4,ARG5,ARG6);				\
  }
#define WARNING7(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7)	\
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7);			\
  }
#define WARNING8(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7,ARG8) \
  {									\
    message2_(stdout,"WARNING",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7,ARG8);			\
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
    terminate_();							\
  }
#define ERROR1(FUNCTION,MESSAGE,ARG1)					\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1);							\
    terminate_();							\
  }
#define ERROR2(FUNCTION,MESSAGE,ARG1,ARG2)				\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2);						\
    terminate_();							\
  }
#define ERROR3(FUNCTION,MESSAGE,ARG1,ARG2,ARG3)				\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3);						\
    terminate_();							\
  }
#define ERROR4(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4)			\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4);					\
    terminate_();							\
  }
#define ERROR5(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5)		\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4,ARG5);				\
    terminate_();							\
  }
#define ERROR6(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6)		\
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4,ARG5,ARG6);				\
    terminate_();							\
  }
#define ERROR8(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7,ARG8) \
  {									\
    message2_(stderr,"ERROR",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
	      ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7,ARG8);			\
    terminate_();							\
  }

//----------------------------------------------------------------------
/// @def      INCOMPLETE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE(FUNCTION)						\
  {									\
    message2_(stdout,"INCOMPLETE",__FILE__,__LINE__,FUNCTION,"");	\
  }

//----------------------------------------------------------------------
/// @def      TRACE
/// @brief    Trace file name and location to stdout
#ifdef CELLO_TRACE

#   define TRACE(MESSAGE)				\
  {							\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",	\
	      MESSAGE);					\
  }

#   define TRACE0					\
  {							\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"", "");	\
  }

#   define TRACE1(MESSAGE,ARG1)				\
  {							\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",	\
	      MESSAGE,ARG1);				\
  }

#   define TRACE2(MESSAGE,ARG1,ARG2)			\
  {							\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",	\
	      MESSAGE,ARG1,ARG2);			\
  }

#   define TRACE3(MESSAGE,ARG1,ARG2,ARG3)		\
  {							\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",	\
	      MESSAGE,ARG1,ARG2,ARG3);			\
  }

#   define TRACE4(MESSAGE,ARG1,ARG2,ARG3,ARG4)		\
  {							\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",	\
	      MESSAGE,ARG1,ARG2,ARG3,ARG4);		\
  }
#   define TRACE5(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5)	\
  {							\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",	\
	      MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5);	\
  }
#   define TRACE6(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6)	\
  {								\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",		\
	      MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6);		\
  }
#   define TRACE7(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7)	\
  {								\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",		\
	      MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7);	\
  }
#   define TRACE8(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7,ARG8)	\
  {								\
    message2_(stdout,"TRACE",__FILE__,__LINE__,"",		\
	      MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7,ARG8);	\
  }
#else

#   define TRACE0				\
  /* This space intentionally left blank */
#   define TRACE(MESSAGE)			\
  /* This space intentionally left blank */
#   define TRACE1(MESSAGE,ARG1)			\
  /* This space intentionally left blank */
#   define TRACE2(MESSAGE,ARG1,ARG2)		\
  /* This space intentionally left blank */
#   define TRACE3(MESSAGE,ARG1,ARG2,ARG3)	\
  /* This space intentionally left blank */
#   define TRACE4(MESSAGE,ARG1,ARG2,ARG3,ARG4)	\
  /* This space intentionally left blank */
#   define TRACE5(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5)	\
  /* This space intentionally left blank */
#   define TRACE6(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6)	\
  /* This space intentionally left blank */
#   define TRACE7(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7)	\
  /* This space intentionally left blank */
#   define TRACE8(MESSAGE,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7,ARG8)	\
  /* This space intentionally left blank */

#endif

//----------------------------------------------------------------------
/// @def      ASSERT
/// @brief    Equivalent to assert()
#define ASSERT(FUNCTION,MESSAGE,ASSERTION)				\
  {									\
    if (!(ASSERTION)) {							\
      message2_(stderr,"ASSERT",__FILE__,__LINE__,FUNCTION,MESSAGE);	\
      terminate_();							\
    }									\
  }
#define ASSERT1(FUNCTION,MESSAGE,ARG1,ASSERTION)			\
  {									\
    if (!(ASSERTION)) {							\
      message2_(stderr,"ASSERT",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
		ARG1);							\
      terminate_();							\
    }									\
  }
#define ASSERT2(FUNCTION,MESSAGE,ARG1,ARG2,ASSERTION)			\
  {									\
    if (!(ASSERTION)) {							\
      message2_(stderr,"ASSERT",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
		ARG1,ARG2);						\
      terminate_();							\
    }									\
  }
#define ASSERT3(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ASSERTION)		\
  {									\
    if (!(ASSERTION)) {							\
      message2_(stderr,"ASSERT",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
		ARG1,ARG2,ARG3);					\
      terminate_();							\
    }									\
  }
#define ASSERT4(FUNCTION,MESSAGE,ARG1,ARG2,ARG3,ARG4,ASSERTION)		\
  {									\
    if (!(ASSERTION)) {							\
      message2_(stderr,"ASSERT",__FILE__,__LINE__,FUNCTION,MESSAGE,	\
		ARG1,ARG2,ARG3,ARG4);					\
      terminate_();							\
    }									\
  }

extern void message2_
(FILE * fp,
 const char * type, 
 const char * file, 
 int line, 
 const char * function, 
 const char * message,
 ...);

void terminate_();

#endif /* ERROR_ERROR_HPP */

