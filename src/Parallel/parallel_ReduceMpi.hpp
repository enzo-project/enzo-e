// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Reduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Mar 11 12:07:45 PST 2011
/// @brief    [\ref Parallel] Declaration and implementation of ReduceMpi

#ifndef PARALLEL_REDUCE_MPI_HPP
#define PARALLEL_REDUCE_MPI_HPP

#include "error.hpp"

class ReduceMpi : public Reduce {

  /// @class    ReduceMpi
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Implementation of ReduceMpi

public: // interface

  /// Constructor
  ReduceMpi(GroupProcess * group_process) throw()
    : Reduce (group_process)
  { /* EMPTY */ };
    
  /// Destructor
  ~ReduceMpi() throw()
  { /* EMPTY */ };

  /// Local reduction of the given value
  virtual int reduce_int
  (  int              local,
     enum_reduce_op   reduce_op )  throw()
  {
#ifdef CONFIG_USE_MPI
    MPI_Datatype mpi_type = MPI_INT;

    MPI_Op mpi_op;
    switch (reduce_op) {
    case (reduce_op_min) :  mpi_op = MPI_MIN; break;
    case (reduce_op_land) : mpi_op = MPI_LAND; break;
    }
    
    GroupProcessMpi * group_mpi = 
      dynamic_cast<GroupProcessMpi *>(group_process_);

    int global;

    MPI_Allreduce (&local, &global, 1, mpi_type, mpi_op, group_mpi->comm());

    return global;
#else
    return 0;
#endif
  }

  /// Local reduction of the given value
  virtual double reduce_double
  (  double              local,
     enum_reduce_op   reduce_op )  throw()
  {
#ifdef CONFIG_USE_MPI
    MPI_Datatype mpi_type = MPI_DOUBLE;

    MPI_Op mpi_op;
    switch (reduce_op) {
    case (reduce_op_min) :  mpi_op = MPI_MIN; break;
    case (reduce_op_land) : mpi_op = MPI_LAND; break;
    }
    
    double global;

    GroupProcessMpi * group_mpi = 
      dynamic_cast<GroupProcessMpi *>(group_process_);

    MPI_Allreduce (&local, &global, 1, mpi_type, mpi_op, group_mpi->comm());

    return global;
#else
    return 0.0;
#endif
  }

};

#endif /* PARALLEL_REDUCE_MPI_HPP */

