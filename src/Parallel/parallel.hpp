// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

/// @file     parallel.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-16
/// @brief    Include file for the Parallel component

enum parallel_enum {
  parallel_serial,
  parallel_mpi
};

#include <string>
#include <vector>
#include <string>
#include <map>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cello.hpp"

#include "error.hpp"

#ifdef CONFIG_USE_MPI
#  include <mpi.h>
#endif

#include "parallel.def"

#include "parallel_Mpi.hpp"

#include "parallel_ParallelGroup.hpp" 

#include "parallel_GroupProcess.hpp"
#ifdef CONFIG_USE_MPI
#   include "parallel_GroupProcessMpi.hpp"
#endif
#include "parallel_GroupProcessSerial.hpp"
#include "parallel_GroupThread.hpp"
#include "parallel_Parallel.hpp"
#include "parallel_Layout.hpp"
#include "parallel_ParallelAffinity.hpp"

#endif /* PARALLEL_HPP */

