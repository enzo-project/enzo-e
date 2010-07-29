// $Id: GroupProcessMpi.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcessMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Jul 23 11:11:07 PDT 2010
/// @brief    Implementation of the GroupProcessMpi class

//----------------------------------------------------------------------

#include "parallel.hpp"

//----------------------------------------------------------------------

GroupProcessMpi::GroupProcessMpi(int process_first,
				 int process_last_plus,
				 int process_stride) throw ()
  : GroupProcess(),
    comm_             (MPI_COMM_WORLD),
    process_first_    (process_first),
    process_last_plus_(process_last_plus),
    process_stride_   (process_stride),
    send_type_        (send_standard)
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ASSERT("GroupProcessMpi::GroupProcessMpi",
	 "process_first out of range",
	 0 <= process_first_ && process_first_ < size);
  ASSERT("GroupProcessMpi::GroupProcessMpi",
	 "process_last_plus out of range",
	 process_last_plus_ == -1 || 
	 ( process_first_ < process_last_plus_ && 
	   process_last_plus_ <= size ));

  // Set default process_last_plus_ to P

  if (process_last_plus_ == -1) process_last_plus_ = size;

  ASSERT("GroupProcessMpi::GroupProcessMpi",
	 "process_stride",
	 1 <= process_stride_ &&
	 process_stride_ <= (process_last_plus_ - process_first_));

  // Compute size_

  size_ = (process_last_plus_ - process_first_) / process_stride_;

  // Update process_last_plus for tight upper-bound

  process_last_plus_ = process_first_ + size_*process_stride_;
  
  rank_ = (rank - process_first_) / process_stride_;

}

//----------------------------------------------------------------------

void GroupProcessMpi::barrier()  throw()
{ 
  MPI_Barrier (comm_); 
};

//----------------------------------------------------------------------

void GroupProcessMpi::sync (int rank, int tag) throw()
{
  char buffer = 1;
  if (rank_ < rank) {
    int ierr = MPI_Send (&buffer, 1, MPI_BYTE, rank, tag, comm_);
    check_mpi_err_("sync",ierr);
  } else if (rank_ > rank) {
    MPI_Status status;
    int ierr = MPI_Recv (&buffer, 1, MPI_BYTE, rank, tag, comm_, &status);
    check_mpi_err_("sync",ierr);
  }
}

//----------------------------------------------------------------------

void * GroupProcessMpi::send_begin 
(int rank_dest, void * buffer, int size, int tag) throw()
{
  int ierr;
  MPI_Request * handle = 0;
  if (send_blocking_) {
    ierr = MPI_Send (buffer, size, MPI_BYTE, rank_dest, tag, comm_);
  } else {
    handle = new MPI_Request;
    ierr = MPI_Isend (buffer, size, MPI_BYTE, rank_dest, tag, comm_,handle);
  }
  check_mpi_err_("send",ierr);
  return (void *)handle;
}

//----------------------------------------------------------------------

bool GroupProcessMpi::send_test (void * handle) throw()
{
  int result = true;
  if (! send_blocking_) {
    MPI_Status status;
    int result;
    int ierr = MPI_Test((MPI_Request*)handle, & result, &status);
    check_mpi_err_("send_test",ierr);
  }
  return result;
}

//----------------------------------------------------------------------

void GroupProcessMpi::send_wait (void * handle) throw()
{
  if (! send_blocking_) {
    MPI_Status status;
    int ierr = MPI_Wait((MPI_Request*)handle, &status);
    check_mpi_err_("send_wait",ierr);
  }
}

//----------------------------------------------------------------------

void GroupProcessMpi::send_end (void * handle) throw()
{
  delete (MPI_Request *) handle;
}

//----------------------------------------------------------------------

void * GroupProcessMpi::recv_begin 
(int rank_source, void * buffer, int size, int tag) throw()
{
  int ierr;
  MPI_Request * handle = 0;
  if (recv_blocking_) {
    MPI_Status status;
    ierr = MPI_Recv (buffer, size, MPI_BYTE, rank_source, tag, comm_, &status);
  } else {
    handle = new MPI_Request;
    ierr = MPI_Irecv (buffer, size, MPI_BYTE, rank_source, tag, comm_, handle);
  }
  check_mpi_err_("recv",ierr);
  return handle;
}

//----------------------------------------------------------------------

bool GroupProcessMpi::recv_test (void * handle) throw()
{
  int result = true;
  if (! recv_blocking_) {
    MPI_Status status;
    int result;
    int ierr = MPI_Test((MPI_Request*)handle, &result, &status);
    check_mpi_err_("recv_test",ierr);
  }
  return result;
}

//----------------------------------------------------------------------

void GroupProcessMpi::recv_wait(void * handle) throw()
{
  if (! recv_blocking_) {
    MPI_Status status;
    int ierr = MPI_Wait((MPI_Request*) handle, &status);
    check_mpi_err_("recv_wait",ierr);
  }
}

//----------------------------------------------------------------------

void GroupProcessMpi::recv_end (void * handle) throw()
{
  delete (MPI_Request *) handle;
}

//----------------------------------------------------------------------

void GroupProcessMpi::bulk_send_add(int rank_dest, void * buffer, int size, int tag) throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessMpi::bulk_send() throw()
{
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessMpi::bulk_send_wait(void * handle) throw()
{
}

//----------------------------------------------------------------------

void GroupProcessMpi::bulk_recv_add(int rank_source, void * buffer, int size, int tag) throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessMpi::bulk_recv() throw()
{
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessMpi::bulk_recv_wait(void * handle) throw()
{
}

//----------------------------------------------------------------------

