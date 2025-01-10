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
// ENUMERATION DECLARATIONS
//----------------------------------------------------------------------

/// @enum     counter_type_enum
/// @brief    Counter value type
enum counter_type_enum {
  counter_type_unknown,
  counter_type_rel,
  counter_type_abs,
  counter_type_papi,
  counter_type_user,
  num_counter_type
};

enum index_enum {
  perf_index_time,
  perf_index_bytes,
  perf_index_bytes_high,
  perf_index_bytes_highest,
  perf_index_bytes_available,
  perf_index_last,
  num_perf_index = perf_index_last
};

/// @enum    perf_region
/// @brief   region ID's for the Simulation performance object
enum perf_region {
  perf_unknown,
  perf_simulation,
  perf_cycle,
  perf_initial,
  perf_adapt,
  perf_adapt_post,
  perf_adapt_enter,
  perf_adapt_enter_post,
  perf_adapt_end,
  perf_adapt_end_post,
  perf_adapt_update,
  perf_adapt_update_post,
  perf_adapt_next,
  perf_adapt_next_post,
  perf_adapt_called,
  perf_adapt_called_post,
  perf_adapt_exit,
  perf_adapt_exit_post,
  perf_adapt_delete,
  perf_adapt_delete_post,
  perf_adapt_recv_level,
  perf_adapt_recv_level_post,
  perf_adapt_recv_child,
  perf_adapt_recv_child_post,
  perf_reduce,
  perf_reduce_stopping,
  perf_reduce_adapt,
  perf_reduce_charm,
  perf_reduce_initialize,
  perf_reduce_output,
  perf_reduce_restart,
  perf_reduce_balance,
  perf_reduce_method_balance,
  perf_reduce_method_check,
  perf_reduce_method_inference,
  perf_reduce_method_m1_closure,
  perf_reduce_method_turbulence,
  perf_reduce_simulation,
  perf_reduce_solver_bicgstab,
  perf_reduce_solver_cg,
  perf_reduce_solver_dd,
  perf_reduce_solver_mg0,
  perf_reduce_method_debug,
  perf_reduce_method_flux_correct,
  perf_reduce_method_order_hilbert,
  perf_reduce_method_order_morton,
  perf_reduce_method_output,
  perf_refresh,
  perf_refresh_post,
  perf_refresh_start,
  perf_refresh_start_post,
  perf_refresh_recv,
  perf_refresh_recv_post,
  perf_refresh_child,
  perf_refresh_child_post,
  perf_refresh_exit,
  perf_refresh_exit_post,
#ifdef CONFIG_SMP_MODE
  perf_smp,
  perf_smp_field_face,
  perf_smp_hierarchy,
  perf_smp_initial_music,
  perf_smp_initial_value,
  perf_smp_method_close_files,
  perf_smp_solver_bcg,
#endif
  perf_balance,
  perf_control,
  perf_method,
  perf_solver,
  perf_output,
  perf_stopping,
  perf_block,
  perf_exit,
#ifdef CONFIG_USE_GRACKLE
  perf_grackle,
#endif
  num_perf_region
};

//----------------------------------------------------------------------
// MACRO DECLARATIONS
//----------------------------------------------------------------------

#ifdef CONFIG_USE_PERFORMANCE
#   define PERF_START(INDEX)                                    \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__)
#   define PERF_STOP(INDEX)                                     \
  cello::performance()->stop_region(INDEX,__FILE__,__LINE__)

#   define PERF_ADAPT_START(INDEX)                                      \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__);          \
  cello::performance()->start_region(perf_adapt,__FILE__,__LINE__);
#   define PERF_ADAPT_STOP(INDEX)                                       \
  cello::performance()->stop_region(INDEX,__FILE__,__LINE__);           \
  cello::performance()->stop_region(perf_adapt,__FILE__,__LINE__)
#   define PERF_ADAPT_POST(INDEX)                                       \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__);          \
  cello::performance()->start_region(perf_adapt_post,__FILE__,__LINE__);

#   define PERF_REDUCE_START(INDEX)                                     \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__);          \
  cello::performance()->start_region(perf_reduce,__FILE__,__LINE__)
#   define PERF_REDUCE_STOP(INDEX)                                      \
  cello::performance()->stop_region(INDEX,__FILE__,__LINE__);           \
  cello::performance()->stop_region(perf_reduce,__FILE__,__LINE__)

#   define PERF_REFRESH_START(INDEX)                                      \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__)
#   define PERF_REFRESH_STOP(INDEX)                                       \
  cello::performance()->stop_region(INDEX,__FILE__,__LINE__)
#   define PERF_REFRESH_POST(INDEX)                                       \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__)

#   define PERF_SOLVER_START(SOLVER)                                    \
  cello::performance()->start_region((SOLVER)->index_perf(),__FILE__,__LINE__); \
  cello::performance()->start_region(perf_solver)
#   define PERF_SOLVER_STOP(SOLVER)                                     \
  cello::performance()->stop_region((SOLVER)->index_perf(),__FILE__,__LINE__); \
  cello::performance()->stop_region(perf_solver)
#   define PERF_METHOD_START(METHOD)                                    \
  cello::performance()->start_region((METHOD)->index_perf(),__FILE__,__LINE__); \
  cello::performance()->start_region(perf_method)
#   define PERF_METHOD_STOP(METHOD)                                     \
  cello::performance()->stop_region((METHOD)->index_perf(),__FILE__,__LINE__); \
  cello::performance()->stop_region(perf_method)
#ifdef CONFIG_SMP_MODE
#   define PERF_SMP_START(INDEX)                                \
  cello::performance()->start_region(INDEX,__FILE__,__LINE__);  \
  cello::performance()->start_region(perf_smp)
#   define PERF_SMP_STOP(INDEX)                                 \
  cello::performance()->stop_region(INDEX,__FILE__,__LINE__);   \
  cello::performance()->stop_region(perf_smp)
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
