// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Mpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Oct 15 10:40:37 PDT 2009 
/// @brief    [\ref Parallel] Interface for the Mpi class

#ifndef PARALLEL_MPI_HPP
#define PARALLEL_MPI_HPP

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
#ifdef CONFIG_USE_MPI
    MPI_Init(argc,argv);
#endif
  }

  /// Finalize MPI
  static void finalize()
  {
#ifdef CONFIG_USE_MPI
    MPI_Finalize();
#endif
  };

  /// Abort execution abruptly
  static void abort()
  {
#ifdef CONFIG_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#else
    exit(1);
#endif
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
#ifdef CONFIG_USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    size = 1;
#endif
    return size;
  };

  /// Get MPI rank
  static int rank()
  { 
    int rank; 
#ifdef CONFIG_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
#else
    rank = 0;
#endif
    return rank; 
  };

  /// Exit the program
  static void barrier()
  {
#ifdef CONFIG_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  /// Return whether this is the root process
  static bool is_root() 
  { 
    return rank() == 0; 
  };

};

#endif /* PARALLEL_MPI_HPP */

