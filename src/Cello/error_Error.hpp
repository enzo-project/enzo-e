// See LICENSE_CELLO file for license and copyright information

/// @file     error_Error.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-07
/// @brief    Declaration of the Error class

#include "charm++.h"

#ifndef ERROR_ERROR_HPP
#define ERROR_ERROR_HPP

//----------------------------------------------------------------------
/// @def      ERROR_LENGTH
/// @brief    Maximum length of error and warning messages
#define ERROR_LENGTH 255

//----------------------------------------------------------------------
/// @def      WARNING
/// @brief    Handle a (non-lethal) warning message
#define WARNING(F,M)                                    \
  { cello::message                                      \
      (stdout,"WARNING",__FILE__,__LINE__,F,M); }
#define WARNING1(F,M,A1)                                \
  { cello::message                                      \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1); }
#define WARNING2(F,M,A1,A2)                                     \
  { cello::message                                              \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2);  }
#define WARNING3(F,M,A1,A2,A3)                                  \
  { cello::message                                              \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3);  }
#define WARNING4(F,M,A1,A2,A3,A4)                               \
  { cello::message                                              \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4);  }
#define WARNING5(F,M,A1,A2,A3,A4,A5)                                    \
  { cello::message                                                      \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5);  }
#define WARNING6(F,M,A1,A2,A3,A4,A5,A6)                                 \
  { cello::message                                                      \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6);  }
#define WARNING7(F,M,A1,A2,A3,A4,A5,A6,A7)                              \
  { cello::message                                                      \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7);  }
#define WARNING8(F,M,A1,A2,A3,A4,A5,A6,A7,A8)                           \
  { cello::message                                                      \
      (stdout,"WARNING",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7,A8); }

//----------------------------------------------------------------------
/// @def      ERROR
/// @brief    Handle a (lethal) error message
#define ERROR(F,M)                              \
  { cello::message                              \
      (stderr,"ERROR",__FILE__,__LINE__,F,M);   \
    cello::error(); }
#define ERROR1(F,M,A1)                                  \
  { cello::message                                      \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1);        \
    cello::error(); }
#define ERROR2(F,M,A1,A2)                               \
  { cello::message                                      \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2);     \
    cello::error(); }
#define ERROR3(F,M,A1,A2,A3)                            \
  { cello::message                                      \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3);  \
    cello::error(); }
#define ERROR4(F,M,A1,A2,A3,A4)                                 \
  { cello::message                                              \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4);       \
    cello::error(); }
#define ERROR5(F,M,A1,A2,A3,A4,A5)                              \
  { cello::message                                              \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5);    \
    cello::error(); }
#define ERROR6(F,M,A1,A2,A3,A4,A5,A6)                           \
  { cello::message                                              \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6); \
    cello::error(); }
#define ERROR7(F,M,A1,A2,A3,A4,A5,A6,A7)                                \
  { cello::message                                                      \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7);      \
    cello::error(); }
#define ERROR8(F,M,A1,A2,A3,A4,A5,A6,A7,A8)				\
  { cello::message                                                      \
      (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7,A8);   \
    cello::error(); }

//----------------------------------------------------------------------
/// @def      INCOMPLETE
/// @brief    Placeholder for code that is incomplete
#define INCOMPLETE(M)                                                   \
  { cello::message(stdout,"INCOMPLETE",__FILE__,__LINE__,"",M);  }

//----------------------------------------------------------------------
/// @def      TRACE
/// @brief    Trace file name and location to stdout

#ifdef CELLO_TRACE

#   define TRACE(M)                                     \
  { cello::message                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M); }
#   define TRACE0                                       \
  { cello::message                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", ""); }
#   define TRACE1(M,A1)                                 \
  { cello::message                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1); }
#   define TRACE2(M,A1,A2)                              \
  { cello::message                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2); }
#   define TRACE3(M,A1,A2,A3)                                   \
  { cello::message                                              \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3); }
#   define TRACE4(M,A1,A2,A3,A4)                                \
  { cello::message                                              \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4); }
#   define TRACE5(M,A1,A2,A3,A4,A5)					\
  { cello::message                                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5); }
#   define TRACE6(M,A1,A2,A3,A4,A5,A6)					\
  { cello::message                                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6); }
#   define TRACE7(M,A1,A2,A3,A4,A5,A6,A7)				\
  { cello::message                                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7); }
#   define TRACE8(M,A1,A2,A3,A4,A5,A6,A7,A8)				\
  { cello::message                                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8); }
#   define TRACE9(M,A1,A2,A3,A4,A5,A6,A7,A8,A9)				\
  { cello::message                                                      \
      (stdout,"TRACE",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8,A9); }

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
#   define TRACE9(M,A1,A2,A3,A4,A5,A6,A7,A8,A9) \
  /* This space intentionally left blank */


#endif /* CELLO_TRACE */

#ifdef CELLO_TRACE_CHARM

#   define TRACEPUP                                             \
  {                                                             \
    CkPrintf("%d TRACEPUP %s:%d %p [%s] [%s]\n",       \
             CkMyPe(),__FILE__,__LINE__,                        \
             (void*)this,                                       \
             p.isPacking()   ? "Packing"   :                    \
             (p.isUnpacking() ? "Unpacking" :                   \
              (p.isSizing()    ? "Sizing"    :"UNKNOWN")),      \
             (p.isDeleting()  ? "Deleting" : "Not deleting"));  \
           fflush(stdout);                                      \
           }


#else
#   define TRACEPUP                             \
  /* This space intentionally left blank */
#endif


//----------------------------------------------------------------------
/// @def      DEBUG
/// @brief    Write debug statements to stderr
#ifdef CELLO_DEBUG

#   define DEBUG(M)                                     \
  {  cello::message                                     \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M); }
#   define DEBUG0                                       \
  { cello::message                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", ""); }
#   define DEBUG1(M,A1)                                 \
  { cello::message                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1); }
#   define DEBUG2(M,A1,A2)                              \
  { cello::message                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2); }
#   define DEBUG3(M,A1,A2,A3)                                   \
  { cello::message                                              \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3); }
#   define DEBUG4(M,A1,A2,A3,A4)                                \
  { cello::message                                              \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4); }
#   define DEBUG5(M,A1,A2,A3,A4,A5)					\
  { cello::message                                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5); }
#   define DEBUG6(M,A1,A2,A3,A4,A5,A6)					\
  { cello::message                                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6); }
#   define DEBUG7(M,A1,A2,A3,A4,A5,A6,A7)				\
  { cello::message                                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7); }
#   define DEBUG8(M,A1,A2,A3,A4,A5,A6,A7,A8)				\
  { cello::message                                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8); }
#   define DEBUG9(M,A1,A2,A3,A4,A5,A6,A7,A8,A9)				\
  { cello::message                                                      \
      (stderr,"DEBUG",__FILE__,__LINE__,"", M,A1,A2,A3,A4,A5,A6,A7,A8,A9); }

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

#define ASSERT(F,M,A)                           \
  {  if (!(A))                                  \
      { cello::message                          \
        (stderr,"ERROR",__FILE__,__LINE__,F,M); \
        cello::error(); } }
#define ASSERT1(F,M,A1,A)                               \
  {  if (!(A))                                          \
      { cello::message                                  \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1);      \
        cello::error(); } }
#define ASSERT2(F,M,A1,A2,A)                            \
  {  if (!(A))                                          \
      { cello::message                                  \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2);   \
        cello::error(); } }
#define ASSERT3(F,M,A1,A2,A3,A)                                 \
  {  if (!(A))                                                  \
      { cello::message                                          \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3);        \
        cello::error(); } }
#define ASSERT4(F,M,A1,A2,A3,A4,A)                              \
  {  if (!(A))                                                  \
      { cello::message                                          \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4);     \
        cello::error(); } }
#define ASSERT5(F,M,A1,A2,A3,A4,A5,A)                           \
  {  if (!(A))                                                  \
      { cello::message                                          \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5);  \
        cello::error(); } }
#define ASSERT6(F,M,A1,A2,A3,A4,A5,A6,A)				\
  {  if (!(A))                                                          \
      { cello::message                                                  \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6);       \
        cello::error(); } }
#define ASSERT7(F,M,A1,A2,A3,A4,A5,A6,A7,A)				\
  {  if (!(A))                                                          \
      { cello::message                                                  \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7);    \
        cello::error(); } }
#define ASSERT8(F,M,A1,A2,A3,A4,A5,A6,A7,A8,A)				\
  {  if (!(A))                                                          \
      { cello::message                                                  \
        (stderr,"ERROR",__FILE__,__LINE__,F,M,A1,A2,A3,A4,A5,A6,A7,A8); \
        cello::error(); } }


namespace cello {
  extern void message
  (FILE * fp,
   const char * type, 
   const char * file, 
   int line, 
   const char * function, 
   const char * message,
   ...);

  void error();
}

#endif /* ERROR_ERROR_HPP */

