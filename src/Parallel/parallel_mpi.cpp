// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      mpi.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Oct 15 10:41:28 PDT 2009
/// @brief     MPI helper functions

#include <mpi.h>

#include "cello.h"

#include "parallel_mpi.hpp"
 
//====================================================================
// PUBLIC FUNCTIONS
//====================================================================

ParallelMpi::ParallelMpi()
  : size_(1),
    rank_(0),
    send_blocking_(true),
    recv_blocking_(true),
    send_type_(send_type_standard)
/**
 *********************************************************************
 *
 * @param         none
 * @return        none
 *
 * Create an ParallelMpi object
 *
 *********************************************************************
 */
{
}

ParallelMpi::~ParallelMpi()
/**
 *********************************************************************
 *
 * @param         none
 * @return        none
 *
 * Delete an ParallelMpi object
 *
 *********************************************************************
 */
{
}
    
//====================================================================
// PRIVATE FUNCTIONS
//====================================================================

