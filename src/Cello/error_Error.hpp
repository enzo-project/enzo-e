// See LICENSE_CELLO file for license and copyright information

/// @file     error_Error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-07
/// @brief    Declaration of the Error class

#ifndef ERROR_ERROR_HPP
#define ERROR_ERROR_HPP

//----------------------------------------------------------------------
/// @def      ERROR_LENGTH
/// @brief    Maximum length of error and warning messages
#define ERROR_LENGTH 255

//----------------------------------------------------------------------
/// @def      WARNING
/// @brief    Handle a (non-lethal) warning message
#define WARNING(F,M)					\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M); }
#define WARNING1(F,M,A1)				\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1); }
#define WARNING2(F,M,A1,A2)					\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2);  }
#define WARNING3(F,M,A1,A2,A3)					\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3);  }
#define WARNING4(F,M,A1,A2,A3,A4)				\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4);  }
#define WARNING5(F,M,A1,A2,A3,A4,A5)					\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5);  }
#define WARNING6(F,M,A1,A2,A3,A4,A5,A6)					\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6);  }
#define WARNING7(F,M,A1,A2,A3,A4,A5,A6,A7)				\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7);  }
#define WARNING8(F,M,A1,A2,A3,A4,A5,A6,A7,A8)				\
  { m2_(stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7,A8); }

//----------------------------------------------------------------------
/// @def      UNTESTED
/// @brief    Handle a (non-lethal) untested message
#define UNTESTED(F)					\
  { m2_(stdout,"UNTESTED",__FILE__,__LINE__,F,""); }

//----------------------------------------------------------------------
/// @def      ERROR
/// @brief    Handle a (lethal) error message
#define ERROR(F,M)					\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M); t_(); }
#define ERROR1(F,M,A1)						\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1); t_(); }
#define ERROR2(F,M,A1,A2)					\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2); t_(); }
#define ERROR3(F,M,A1,A2,A3)					\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3); t_(); }
#define ERROR4(F,M,A1,A2,A3,A4)						\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4); t_(); }
#define ERROR5(F,M,A1,A2,A3,A4,A5)					\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5); t_(); }
#define ERROR6(F,M,A1,A2,A3,A4,A5,A6)					\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6); t_(); }
#define ERROR8(F,M,A1,A2,A3,A4,A5,A6,A7,A8)				\
  { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7,A8); t_(); }

//----------------------------------------------------------------------
/// @def      INCOMPLETE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE(M)					\
  { m2_(stdout,"INCOMPLETE",__FILE__,__LINE__,"",M);  }
#define INCOMPLETE1(M,A1)					\
  { m2_(stdout,"INCOMPLETE",__FILE__,__LINE__,"",M,A1);  }

//----------------------------------------------------------------------
/// @def      TRACE
/// @brief    Trace file name and location to stdout
#ifdef CELLO_TRACE

#   define TRACE(M)					\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M); }
#   define TRACE0						\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", ""); }
#   define TRACE1(M,A1)					\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1); }
#   define TRACE2(M,A1,A2)					\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2); }
#   define TRACE3(M,A1,A2,A3)					\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3); }
#   define TRACE4(M,A1,A2,A3,A4)				\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4); }
#   define TRACE5(M,A1,A2,A3,A4,A5)					\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5); }
#   define TRACE6(M,A1,A2,A3,A4,A5,A6)					\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6); }
#   define TRACE7(M,A1,A2,A3,A4,A5,A6,A7)				\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7); }
#   define TRACE8(M,A1,A2,A3,A4,A5,A6,A7,A8)				\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8); }
#   define TRACE9(M,A1,A2,A3,A4,A5,A6,A7,A8,A9)				\
  { m2_(stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8,A9); }
#   define TRACEPUP							\
  { m2_(stdout,"TRACEPUP",__FILE__,__LINE__,"",				\
	p.isPacking()?"isPacking":(p.isUnpacking()?"isUnpacking":"isSizing")); }

#else /* CELLO_TRACE */

#   define TRACE0				\
  /* This space intentionally left blank */
#   define TRACE(M)				\
  /* This space intentionally left blank */
#   define TRACE1(M,A1)				\
  /* This space intentionally left blank */
#   define TRACE2(M,A1,A2)			\
  /* This space intentionally left blank */
#   define TRACE3(M,A1,A2,A3)			\
  /* This space intentionally left blank */
#   define TRACE4(M,A1,A2,A3,A4)		\
  /* This space intentionally left blank */
#   define TRACE5(M,A1,A2,A3,A4,A5)		\
  /* This space intentionally left blank */
#   define TRACE6(M,A1,A2,A3,A4,A5,A6)		\
  /* This space intentionally left blank */
#   define TRACE7(M,A1,A2,A3,A4,A5,A6,A7)	\
  /* This space intentionally left blank */
#   define TRACE8(M,A1,A2,A3,A4,A5,A6,A7,A8)	\
  /* This space intentionally left blank */
#   define TRACE9(M,A1,A2,A3,A4,A5,A6,A7,A8,A9)				\
  /* This space intentionally left blank */

#   define TRACEPUP							\
  /* This space intentionally left blank */

#endif /* CELLO_TRACE */

#ifdef CELLO_TRACE_CHARM
#   define TRACE_CHARM(M)				\
  { m2_(stdout,"TRACE_CHARM",__FILE__,__LINE__,"",M); }
#else
#   define TRACE_CHARM(M)				\
  /* This space intentionally left blank */
#endif


//----------------------------------------------------------------------
/// @def      DEBUG
/// @brief    Write debug statements to stderr
#ifdef CELLO_DEBUG

#   define DEBUG(M)					\
  {  m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M); }
#   define DEBUG0						\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", ""); }
#   define DEBUG1(M,A1)					\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1); }
#   define DEBUG2(M,A1,A2)					\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2); }
#   define DEBUG3(M,A1,A2,A3)					\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3); }
#   define DEBUG4(M,A1,A2,A3,A4)				\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4); }
#   define DEBUG5(M,A1,A2,A3,A4,A5)					\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5); }
#   define DEBUG6(M,A1,A2,A3,A4,A5,A6)					\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6); }
#   define DEBUG7(M,A1,A2,A3,A4,A5,A6,A7)				\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7); }
#   define DEBUG8(M,A1,A2,A3,A4,A5,A6,A7,A8)				\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8); }
#   define DEBUG9(M,A1,A2,A3,A4,A5,A6,A7,A8,A9)				\
  { m2_(stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8,A9); }

#else /* CELLO_DEBUG */

#   define DEBUG0				\
  /* This space intentionally left blank */
#   define DEBUG(M)				\
  /* This space intentionally left blank */
#   define DEBUG1(M,A1)				\
  /* This space intentionally left blank */
#   define DEBUG2(M,A1,A2)			\
  /* This space intentionally left blank */
#   define DEBUG3(M,A1,A2,A3)			\
  /* This space intentionally left blank */
#   define DEBUG4(M,A1,A2,A3,A4)		\
  /* This space intentionally left blank */
#   define DEBUG5(M,A1,A2,A3,A4,A5)		\
  /* This space intentionally left blank */
#   define DEBUG6(M,A1,A2,A3,A4,A5,A6)		\
  /* This space intentionally left blank */
#   define DEBUG7(M,A1,A2,A3,A4,A5,A6,A7)	\
  /* This space intentionally left blank */
#   define DEBUG8(M,A1,A2,A3,A4,A5,A6,A7,A8)	\
  /* This space intentionally left blank */
#   define DEBUG9(M,A1,A2,A3,A4,A5,A6,A7,A8,A9)	\
  /* This space intentionally left blank */

#endif /* CELLO_DEBUG */

//----------------------------------------------------------------------
/// @def      ASSERT
/// @brief    Equivalent to assert()

// #ifdef CELLO_DEBUG

#define ASSERT(F,M,A)							\
  {  if (!(A)) { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M); t_(); } }
#define ASSERT1(F,M,A1,A)						\
  {  if (!(A)) { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1); t_(); } }
#define ASSERT2(F,M,A1,A2,A)						\
  {  if (!(A)) { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2); t_(); } }
#define ASSERT3(F,M,A1,A2,A3,A)						\
  {  if (!(A)) { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3); t_(); } }
#define ASSERT4(F,M,A1,A2,A3,A4,A)					\
  {  if (!(A)) { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4); t_(); } }
#define ASSERT5(F,M,A1,A2,A3,A4,A5,A)					\
  {  if (!(A)) { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5); t_(); } }
#define ASSERT6(F,M,A1,A2,A3,A4,A5,A6,A)				\
  {  if (!(A)) { m2_(stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6); t_(); } }

// #else  /* CELLO_DEBUG */

// #define ASSERT(F,M,A) /* NULL */
// #define ASSERT1(F,M,A1,A)  /* NULL */
// #define ASSERT2(F,M,A1,A2,A) /* NULL */
// #define ASSERT3(F,M,A1,A2,A3,A) /* NULL */
// #define ASSERT4(F,M,A1,A2,A3,A4,A) /* NULL */
// #define ASSERT5(F,M,A1,A2,A3,A4,A5,A) /* NULL */

// #endif /* CELLO_DEBUG */

extern void m2_
(FILE * fp,
 const char * type, 
 const char * file, 
 int line, 
 const char * function, 
 const char * message,
 ...);

void t_();

#endif /* ERROR_ERROR_HPP */

