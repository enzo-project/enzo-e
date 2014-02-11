// See LICENSE_CELLO file for license and copyright information

/// @file     _comm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-27
/// @brief    Private include file for the \ref Comm component 

#ifndef _COMM_HPP
#define _COMM_HPP

/// @enum     phase_sync_type
/// @brief    adapt phase for p_join()
enum phase_sync_type {
  phase_sync_unknown,
  phase_sync_adapt_called,
  phase_sync_adapt_next,
  phase_sync_adapt_end,
  phase_sync_refresh
};
#define PHASE_SYNC_SIZE 5

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include "pup_stl.h"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "comm_CommBlock.hpp"

#endif /* _COMM_HPP */
