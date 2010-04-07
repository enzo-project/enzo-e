// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      mpi.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Oct 15 10:41:28 PDT 2009
/// @brief     MPI helper functions

#include <mpi.h>

#include "cello.h"

#include "parallel.hpp"

ParallelMpi ParallelMpi::instance_; // (singleton design pattern)
