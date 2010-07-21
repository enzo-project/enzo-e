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

#include "parallel_Parallel.hpp"
#include "parallel_ParallelMpi.hpp"
#include "parallel_ParallelSerial.hpp"
#include "parallel_ParallelCreate.hpp"

#include "parallel_ParallelGroup.hpp"
#include "parallel_ParallelLayout.hpp"
#include "parallel_ParallelAffinity.hpp"

#endif /* PARALLEL_HPP */

