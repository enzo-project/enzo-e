// $Id: parallel_affinity.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef AFFINITY_HPP
#define AFFINITY_HPP

/// @file     parallel_affinity.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 16:40:35 PDT 2010
/// @brief    Declaration of the Affinity class

class Affinity {

  /// @class    Affinity
  /// @ingroup  Parallel
  /// @brief    Generalization of MPI rank for MPI + OMP + UPC

public: // interface

  /// Initialize the Affinity object
  Affinity();

  /// Delete the Affinity object
  ~Affinity();

private: // functions


private: // attributes

  int mpi_size_;
  int mpi_rank_;
  int omp_size_;
  int omp_rank_;
  int upc_size_;
  int upc_rank_;

};

#endif /* AFFINITY_HPP */

