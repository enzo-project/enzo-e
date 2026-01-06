// See LICENSE_CELLO file for license and copyright information

/// @file     _performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    Private include file for the \ref Performance component

#ifndef _PERFORMANCE_HPP
#define _PERFORMANCE_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <vector>
#include <map>
#include <stack>
#include <string>
#include <sstream>
#include <sys/resource.h>

#ifdef __linux__
#   include <unistd.h>
#endif
#ifdef CONFIG_USE_PAPI
#  include "papi.h"
#endif

//----------------------------------------------------------------------
// MACRO DECLARATIONS
//----------------------------------------------------------------------

#ifdef CONFIG_USE_PERFORMANCE
#   define PERF_START(INDEX)                                    \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__)
#   define PERF_STOP(INDEX)                                     \
  cello::performance()->stop_region(INDEX,__FILE__,__LINE__)

#   define PERF_ADAPT_START(INDEX)              \
  PERF_START(INDEX);                            \
  PERF_START(perf_rindex_adapt);
#   define PERF_ADAPT_STOP(INDEX)               \
  PERF_STOP(perf_rindex_adapt);                 \
  PERF_STOP(INDEX);
#   define PERF_ADAPT_POST(INDEX)               \
  PERF_START(INDEX);                            \
  PERF_START(perf_rindex_adapt_post);

#   define PERF_REDUCE_START(INDEX)             \
  PERF_START(perf_rindex_reduce);               \
  PERF_START(INDEX);
#   define PERF_REDUCE_STOP(INDEX)              \
  PERF_STOP(INDEX);                             \
  PERF_STOP(perf_rindex_reduce);

#   define PERF_REFRESH_START(INDEX)            \
  PERF_START(INDEX);
#   define PERF_REFRESH_STOP(INDEX)             \
  PERF_STOP(INDEX);
#   define PERF_REFRESH_POST(INDEX)             \
  PERF_START(INDEX);

#   define PERF_SOLVER_START(SOLVER)            \
  PERF_START((SOLVER)->index_perf());           \
  PERF_START(perf_rindex_solver);
#   define PERF_SOLVER_STOP(SOLVER)             \
  PERF_STOP(perf_rindex_solver);                \
  PERF_STOP((SOLVER)->index_perf());

#   define PERF_METHOD_START(METHOD)            \
  PERF_START((METHOD)->index_perf());           \
  PERF_START(perf_rindex_method);
#   define PERF_METHOD_STOP(METHOD)             \
  PERF_STOP(perf_rindex_method);                \
  PERF_STOP((METHOD)->index_perf());

#ifdef CONFIG_SMP_MODE
#   define PERF_SMP_START(INDEX)                \
  PERF_START(INDEX);                            \
  PERF_START(perf_rindex_smp);
#   define PERF_SMP_STOP(INDEX)                 \
  PERF_STOP(perf_rindex_smp);                   \
  PERF_STOP(INDEX);
#endif

#else
#   define PERF_START(INDEX) /* ... */
#   define PERF_STOP(INDEX) /* ... */
#   define PERF_ADAPT_START(INDEX) /* ... */
#   define PERF_ADAPT_STOP(INDEX) /* ... */
#   define PERF_ADAPT_POST(INDEX) /* ... */
#   define PERF_REFRESH_START(INDEX) /* ... */
#   define PERF_REFRESH_STOP(INDEX)  /* ... */
#   define PERF_REFRESH_POST(INDEX)  /* ... */
#   define PERF_SOLVER_START(SOLVER) /* ... */
#   define PERF_SOLVER_STOP(SOLVER) /* ... */
#   define PERF_METHOD_START(METHOD) /* ... */
#   define PERF_METHOD_STOP(METHOD) /* ... */
#   define PERF_SMP_START(INDEX) /* ... */
#   define PERF_SMP_STOP(INDEX) /* ... */
#endif

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "performance_Timer.hpp"
#ifdef CONFIG_USE_PAPI  
#include "performance_Papi.hpp"
#endif
#include "performance_Performance.hpp"


#endif /* _PERFORMANCE_HPP */
