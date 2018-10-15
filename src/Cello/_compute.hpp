// See LICENSE_CELLO file for license and copyright information

/// @file     _compute.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Private include file for the \ref Compute component

#ifndef _COMPUTE_HPP
#define _COMPUTE_HPP

//----------------------------------------------------------------------

enum solve_enum {
  solve_unknown,
  solve_leaves,  // Solve on leaf Blocks (default)
  solve_level,   // Solve within a level (e.g. for multigrid smoothers)
  solve_tree,    // Solve in a root-level octree (e.g. for domain decomposition)
  solve_block    // Solve in a block (e.g. MG coarse solve on a single block)
};

extern const char * solve_string[];

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "compute_Compute.hpp"
#include "compute_Matrix.hpp"
#include "compute_Solver.hpp"
#include "compute_SolverNull.hpp"

#endif /* _FIELD_HPP */

