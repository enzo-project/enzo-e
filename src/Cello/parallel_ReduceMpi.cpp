// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_ReduceMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-11-01
/// @brief    Implementation of the ReduceMpi class for parallel reductions
///
/// 

#include "parallel.hpp"

#ifdef CONFIG_USE_MPI

//----------------------------------------------------------------------

int ReduceMpi::reduce_int
(  int              local,
   enum_reduce_op   reduce_op )  throw()
{
#ifdef CONFIG_USE_MPI
  MPI_Datatype mpi_type = MPI_INT;

  MPI_Op mpi_op = get_reduce_type_(reduce_op);

  ASSERT1("ReduceMpi::reduce_int",
	  "Unrecognized operation %d",
	  int(reduce_op),
	  mpi_op != MPI_OP_NULL);

  const GroupProcessMpi * group_mpi = 
    dynamic_cast<const GroupProcessMpi *>(group_process_);

  int global;

  MPI_Allreduce (&local, &global, 1, mpi_type, mpi_op, group_mpi->comm());

  return global;
#else
  return local;
#endif
}

//----------------------------------------------------------------------

double ReduceMpi::reduce_double
(  double              local,
   enum_reduce_op   reduce_op )  throw()
{
#ifdef CONFIG_USE_MPI
  MPI_Datatype mpi_type = MPI_DOUBLE;

  MPI_Op mpi_op = get_reduce_type_(reduce_op);
    
  ASSERT1("ReduceMpi::reduce_int",
	  "Unrecognized operation %d",
	  int(reduce_op),
	  mpi_op != MPI_OP_NULL);

  double global;

  const GroupProcessMpi * group_mpi = 
    dynamic_cast<const GroupProcessMpi *>(group_process_);

  MPI_Allreduce (&local, &global, 1, mpi_type, mpi_op, group_mpi->comm());

  return global;
#else
  return local;
#endif
}

//======================================================================

MPI_Op ReduceMpi::get_reduce_type_(enum_reduce_op reduce_op) const throw()
{
  MPI_Op mpi_op;
  switch (reduce_op) {
  case (reduce_op_min) :  mpi_op = MPI_MIN; break;
  case (reduce_op_max) :  mpi_op = MPI_MAX; break;
  case (reduce_op_sum) :  mpi_op = MPI_SUM; break;
  case (reduce_op_land) : mpi_op = MPI_LAND; break;
  default:                mpi_op = MPI_OP_NULL; break;
  }
  return mpi_op;
}    

#endif /* CONFIG_USE_MPI */
