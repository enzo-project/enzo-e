// See LICENSE_CELLO file for license and copyright information

/// @file     _parallel.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-16
/// @brief    Private include file for the \ref Parallel component

#ifndef _PARALLEL_HPP
#define _PARALLEL_HPP

//----------------------------------------------------------------------
// Enumerations
//----------------------------------------------------------------------

/// @enum enum_reduce_op
/// @bried type of reduction operation
enum enum_reduce_op {
  reduce_op_unknown,
  reduce_op_min,
  reduce_op_land
};

/// @enum enum_reduce_type
/// @bried data type to reduce
enum enum_reduce_type {
  reduce_type_int,
  reduce_type_double
};

enum parallel_enum {
  parallel_serial,
  parallel_mpi
};

#define PROCESS_NULL -1

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>
#include <string>
#include <map>
#include <limits>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef CONFIG_USE_MPI
#  include <mpi.h>
#endif
#ifdef CONFIG_USE_CHARM
#  include "charm++.h"
#endif

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "parallel.def"
#include "parallel_Mpi.hpp"


#include "parallel_GroupProcess.hpp"
#include "parallel_GroupProcessMpi.hpp"
#include "parallel_GroupProcessCharm.hpp"
#include "parallel_GroupProcessSerial.hpp"

//#include "parallel_GroupThread.hpp"

#include "parallel_Reduce.hpp"
#include "parallel_ReduceSerial.hpp"
#include "parallel_ReduceMpi.hpp"
#include "parallel_ReduceCharm.hpp"

// #include "parallel_Parallel.hpp"

#include "parallel_Layout.hpp"


//#include "parallel_ParallelAffinity.hpp"

#endif /* _PARALLEL_HPP */

