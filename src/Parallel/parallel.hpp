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

#include "parallel_Parallel.hpp"
#include "parallel_ParallelMpi.hpp"
#include "parallel_ParallelSerial.hpp"
#include "parallel_ParallelCreate.hpp"

#include "parallel_Layout.hpp"
#include "parallel_Affinity.hpp"

#endif /* PARALLEL_HPP */

