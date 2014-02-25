// See LICENSE_CELLO file for license and copyright information

/// @file     _comm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-27
/// @brief    Private include file for the \ref Comm component 

#ifndef _COMM_HPP
#define _COMM_HPP

//----------------------------------------------------------------------
// Defines
//----------------------------------------------------------------------

// index for child blocks
#define IC3(ic3)  ( ((ic3[0]+2)%2) + 2*( ((ic3[1]+2)%2) + 2*( ((ic3[2]+2)%2) )))

// number of children
#define NC(rank) (1<<(rank))

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

enum phase_type {
  phase_unknown,
  phase_initial,
  phase_adapt,
  phase_compute,
  phase_refresh,
  phase_stopping,
  phase_output
};

enum array_type {
  op_array_unknown,
  op_array_copy,
  op_array_restrict,
  op_array_prolong
};

/// @enum     phase_sync_type
/// @brief    adapt phase for p_join()
enum phase_sync_type {
  phase_sync_unknown,
  phase_sync_adapt_called,
  phase_sync_adapt_enter,
  phase_sync_adapt_next,
  phase_sync_adapt_exit,
  phase_sync_refresh_enter,
  phase_sync_refresh_exit
};
#define PHASE_SYNC_SIZE 7

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <limits>

#include "parallel.def"
#include "mesh.decl.h"
#include "pup_stl.h"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "comm_CommBlock.hpp"

#endif /* _COMM_HPP */
