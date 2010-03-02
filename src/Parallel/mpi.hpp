// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MPI_HPP
#define MPI_HPP

/// @file     mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    Interface for the Mpi class

class Mpi {

  /// @class    Mpi
  /// @ingroup  Parallel
  /// @brief    MPI helper functions

public: // interface

  /// Initialize an Mpi object
  Mpi();

  /// Delete an Mpi object
  ~Mpi();

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

#endif /* MPI_HPP */

