// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    Interface for the Mpi class

#ifndef PARALLEL_MPI_HPP
#define PARALLEL_MPI_HPP

class Mpi {

  /// @class    Mpi
  /// @ingroup  Parallel
  /// @brief    MPI helper functions

private:

  /// Initialize an Mpi object [private to avoid creating a
  /// Mpi object: all functions are static]
  Mpi();

public: // interface

  /// Initialize MPI (virtual)
  static void initialize(int * argc, char ***argv)
  {
    MPI_Init(argc,argv);
  }

  /// Finalize MPI
  static void finalize()
  {
    MPI_Finalize();
  };

  /// Abort execution abruptly
  static void abort()
  {
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /// Exit the program
  static void halt()
  {
    finalize();
    exit (0);
  }

  /// Get MPI size
  static int size()
  {
    int size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  };

  /// Get MPI rank
  static int rank()
  { 
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    return rank; 
  };

  /// Return whether this is the root process
  static bool is_root() 
  { 
    return rank() == 0; 
  };

};

#endif /* PARALLEL_MPI_HPP */

