// See LICENSE_CELLO file for license and copyright information

/// @file     _problem.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Private include file for the \ref Problem component

#ifndef _PROBLEM_HPP
#define _PROBLEM_HPP

/// @enum     sync_enum
/// @brief    synchronization type for refresh
enum sync_enum {
  sync_unknown,   // Unknown synchronization
  sync_none,
  sync_barrier,
  sync_quiescence,
  sync_neighbor,
  sync_face
};

/// @enum     neighbor_enum
/// @brief    neighbor block type
enum neighbor_enum {
  neighbor_unknown, // Unknown neighbor type
  neighbor_leaf,    // neighbors that are leaves, maybe different level
  neighbor_level,   // neighbors is in same level, maybe not leaves
  neighbor_tree     // neighbors that are leaves, but only if in same octree
};
  
//----------------------------------------------------------------------

extern void method_close_files_mutex_init();

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>
#include <limits>
#include <algorithm>

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "problem_Refresh.hpp"
#include "problem_Mask.hpp"
#include "problem_MaskExpr.hpp"
#include "problem_MaskPng.hpp"
#include "problem_ScalarExpr.hpp"
#include "problem_Value.hpp"
//#include "problem_MaskPng.hpp"
//#include "problem_ExprValue.hpp"
#include "problem_Problem.hpp"
#include "problem_Stopping.hpp"
#include "problem_Initial.hpp"
#include "problem_InitialTrace.hpp"
#include "problem_InitialValue.hpp"
#include "problem_Boundary.hpp"
#include "problem_BoundaryPeriodic.hpp"
#include "problem_BoundaryValue.hpp"
#include "problem_Method.hpp"
#include "problem_MethodCloseFiles.hpp"
#include "problem_MethodDebug.hpp"
#include "problem_MethodFluxCorrect.hpp"
#include "problem_MethodNull.hpp"
#include "problem_MethodOrderMorton.hpp"
#include "problem_MethodOutput.hpp"
#include "problem_MethodRefresh.hpp"
#include "problem_MethodTrace.hpp"
#include "problem_Physics.hpp"
#include "problem_Prolong.hpp"
#include "problem_ProlongInject.hpp"
#include "problem_ProlongLinear.hpp"
#include "problem_Restrict.hpp"
#include "problem_RestrictLinear.hpp"
#include "problem_Units.hpp"

#endif /* _PROBLEM_HPP */
