// See LICENSE_CELLO file for license and copyright information

/// @file     _mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Private include file for the \ref Mesh component 

#ifndef _MESH_HPP
#define _MESH_HPP

//----------------------------------------------------------------------
// Defines
//----------------------------------------------------------------------

// index for child blocks
#define IC3(ic3)  ( ((ic3[0]+2)%2) + 2*( ((ic3[1]+2)%2) + 2*( ((ic3[2]+2)%2) )))

// number of children
#define NUM_CHILDREN(rank) (1<<(rank))

// index for neighbors (axis,face)
#define IN(axis,face)  ((face) + 2*(axis))

// index for face (ix,iy,iz)
#define IF3(if3)  ((if3[0]+1) + 3*((if3[1]+1) + 3*(if3[2]+1)))

// index for child ic3[] face if3[]
#define ICF3(ic3,if3)  (IF3(if3) + 27*IC3(ic3))

// number of neighbors
#define NN(rank) (2*(rank))

//----------------------------------------------------------------------
// Enums
//----------------------------------------------------------------------

// KEEP CONSISTENT WITH comm_Block.cpp: phase_string
enum phase_type {
  phase_unknown,
  phase_initial,
  phase_initial_enter = phase_initial,
  phase_initial_exit,
  phase_adapt,
  phase_adapt_enter = phase_adapt,
  phase_adapt_called,
  phase_adapt_next,
  phase_adapt_end,
  phase_adapt_exit,
  phase_compute,
  phase_compute_enter = phase_compute,
  phase_compute_continue,
  phase_compute_exit,
  phase_refresh,
  phase_refresh_enter = phase_refresh,
  phase_refresh_exit,
  phase_stopping,
  phase_stopping_enter = phase_stopping,
  phase_stopping_exit,
  phase_output,
  phase_output_enter = phase_output,
  phase_output_exit,
  phase_restart,
  phase_balance,
  phase_exit,
  phase_last
};

/// @enum   RefreshState
///
/// @brief  New refresh states:
///         RefreshState::inactive (done() called),
///         RefreshState::active   (start() called),
///         RefreshState::ready    (wait() called)

// #define PHASE_COUNT (phase_exit + 1)
// #define PHASE_COUNT 100
extern const char * phase_name[];

/// @enum     adapt_type
/// @brief    Mesh adaptation type: refine, coarsen, or stay the same

enum adapt_type {
  adapt_unknown,
  adapt_coarsen,
  adapt_same,
  adapt_refine
};

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

// Hierarchy and components

class Tree;

#include "mesh_Adapt.hpp"
#include "mesh_Box.hpp"
#include "mesh_Index.hpp"

#include "mesh_Block.hpp"
#include "mesh_Hierarchy.hpp"
#include "mesh_Factory.hpp"

// Tree and components (not used in Cello)
#include "mesh_Node.hpp"
#include "mesh_NodeTrace.hpp"
#include "mesh_BlockTrace.hpp"
#include "mesh_Tree.hpp"
#include "mesh_ItNode.hpp"

// Iterators
#include "mesh_ItFace.hpp"
#include "mesh_ItNeighbor.hpp"
#include "mesh_ItChild.hpp"

// Refinement
#include "mesh_Refine.hpp"
#include "mesh_RefineDensity.hpp"
#include "mesh_RefineMask.hpp"
#include "mesh_RefineShear.hpp"
#include "mesh_RefineSlope.hpp"
#include "mesh_RefineParticleCount.hpp"

#endif /* _MESH_HPP */

