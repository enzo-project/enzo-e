// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    [\ref Parallel] Interface for the Mpi class

#ifndef PARALLEL_MPI_HPP
#define PARALLEL_MPI_HPP

#ifdef CONFIG_USE_MPI

class Mpi {

  /// @class    Mpi
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] MPI helper functions

private:

  /// Initialize an Mpi object [private to avoid creating a
  /// Mpi object: all functions are static]
  Mpi();

public: // interface

  /// Initialize MPI (virtual)
  static void init(int * argc, char ***argv)
  {
    MPI_Init(argc,argv);
  }

  /// Finalize MPI
  static void finalize()
  {
    MPI_Finalize();
  };

  /// Abort execution abruptly
  static void abort(MPI_Comm comm = MPI_COMM_CELLO)
  {
    MPI_Abort(comm,1);
  }

  /// Exit the program
  static void halt()
  {
    finalize();
    exit (0);
  }

  /// Get MPI size
  static int size(MPI_Comm comm = MPI_COMM_CELLO)
  {
    int size; 
    MPI_Comm_size(comm, &size);
    return size;
  };

  /// Get MPI rank
  static int rank(MPI_Comm comm = MPI_COMM_CELLO)
  { 
    int rank; 
    MPI_Comm_rank(comm, &rank); 
    return rank; 
  };

  /// Exit the program
  static void barrier(MPI_Comm comm = MPI_COMM_CELLO)
  {
    MPI_Barrier(comm);
  }

  /// Return whether this is the root process
  static bool is_root() 
  { 
    return rank() == 0; 
  };

};

#endif /* CONFIG_USE_MPI */

#endif /* PARALLEL_MPI_HPP */

