// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

/// @file     parallel.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-16
/// @brief    Include file for the Parallel component

enum parallel_type {
  parallel_serial,
  parallel_mpi
};

#include <string>
#include <vector>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "error.hpp"

#include "parallel_Mpi.hpp"

#include "parallel_Group.hpp"
#include "parallel_GroupProcess.hpp"
#include "parallel_GroupProcessMpi.hpp"
#include "parallel_Layout.hpp"
#include "parallel_Affinity.hpp"

#endif /* PARALLEL_HPP */

