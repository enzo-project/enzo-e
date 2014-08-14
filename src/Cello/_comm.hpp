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
  phase_exit
};

#define PHASE_COUNT (phase_exit + 1)
extern const char * phase_name[];

enum array_type {
  op_array_unknown,
  op_array_copy,
  op_array_restrict,
  op_array_prolong
};

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
