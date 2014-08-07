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

// KEEP CONSISTENT WITH comm_CommBlock.cpp: phase_string
enum phase_type {
  phase_unknown,
  phase_initial,
  phase_adapt,
  phase_compute,
  phase_refresh,
  phase_stopping,
  phase_output,
  phase_restart,
  phase_balance
};

extern const char * phase_name[9];

enum array_type {
  op_array_unknown,
  op_array_copy,
  op_array_restrict,
  op_array_prolong
};

/// @enum     sync_type
/// @brief    adapt phase for p_join()
enum sync_type {
  sync_unknown,
  sync_adapt_enter,
  sync_adapt_called,
  sync_adapt_next,
  sync_adapt_end,
  sync_adapt_exit,
  sync_compute_enter,
  sync_compute_exit,
  sync_output_enter,
  sync_output_exit,
  sync_refresh_enter,
  sync_refresh_exit,
  sync_stopping_enter,
  sync_stopping_exit,
  sync_exit,
};
#define SYNC_SIZE 15

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
