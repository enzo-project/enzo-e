// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_MPI_HPP
#define PARALLEL_MPI_HPP

/// @file     mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    Interface for the ParallelMpi class

class ParallelMpi {

  /// @class    ParallelMpi
  /// @ingroup  Parallel
  /// @brief    MPI helper functions

public: // interface

  /// Initialize an ParallelMpi object
  ParallelMpi();

  /// Delete an ParallelMpi object
  ~ParallelMpi();

private: // attributes

  /// Whether to use blocking sends and receives
  bool blocking_;
  
  /// Whether to use standard, buffered, synchronous, or ready sends
  enum type_send {
    type_standard,
    type_buffered,
    type_synchronous,
    type_ready } type_;

};

#endif /* PARALLEL_MPI_HPP */

