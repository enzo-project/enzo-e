// $Id: parallel_Reduce.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_REDUCE_HPP
#define PARALLEL_REDUCE_HPP

/// @file     parallel_Reduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 10 12:20:48 PST 2011
/// @brief    [\ref Parallel] Declaration and implementation of Reduce

enum enum_mpi_op {
  op_min,
  op_land
};

#include "error.hpp"

template <class T>
class Reduce {

  /// @class    Reduce
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Implementation of reduction operations

public: // interface

  /// Constructor
#ifdef CONFIG_USE_MPI
  Reduce(enum_mpi_op op,MPI_Comm mpi_comm) throw()
    : op_(op),
      mpi_comm_(mpi_comm)
  { 
    reset();
    switch (op_) {
    case op_min:  mpi_op_ = MPI_MIN; break;
    case op_land: mpi_op_ = MPI_LAND; break;
    default:
      ERROR("Reduce::Reduce","Operation not implemented");
    }
  };
#else
  Reduce(enum_mpi_op op) throw()
    : op_(op)
  {
    reset();
  };
#endif
    
  /// Destructor
  ~Reduce() throw()
  {};

  /// Reset to prepare for another reduction
  inline void reset () throw()
  { 
    switch (op_) {
    case op_min: local_ = std::numeric_limits<T>::max(); break;
    case op_land:local_ = true; break;
    default:
      ERROR("Reduce::reset","Operation not implemented");
    }
  };

  /// Local reduction of the given value
  inline void accum (T local) 
  {
    switch (op_) {
    case op_min:  local_ = MIN(local_, local); break;
    case op_land: local_ = local_ && local; break;
    default:
      ERROR("Reduce::accum","Operation not implemented");
    }
  }

  /// Parallel reduction of the stored local value
  inline T reduce () 
  {
    T global;
#ifdef CONFIG_USE_MPI
    MPI_Allreduce (&local_, &global, sizeof(T), MPI_CHAR,mpi_op_,mpi_comm_);
#else
    global = local_;
#endif    
    return global;
  }

private: // attributes

  T local_;
  enum enum_mpi_op op_;
  
#ifdef CONFIG_USE_MPI
  MPI_Comm mpi_comm_;
  MPI_Op mpi_op_;
#endif

private: // attributes

};

#endif /* PARALLEL_REDUCE_HPP */

