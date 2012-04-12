// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcessMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Jul 23 11:11:07 PDT 2010
/// @brief    Implementation of the GroupProcessMpi class


#ifdef CONFIG_USE_MPI

//----------------------------------------------------------------------

#include "cello.hpp"

#include "parallel.hpp"

//----------------------------------------------------------------------

GroupProcessMpi::GroupProcessMpi(int process_first,
				 int process_last_plus) throw ()
  : GroupProcess(),
    comm_             (MPI_COMM_WORLD),
    process_first_    (process_first),
    process_last_plus_(process_last_plus),
    send_type_        (send_standard),
    send_blocking_    (true),
    recv_blocking_    (true)
{
  int size = Mpi::size();
  process_rank_ = Mpi::rank();

  // Check range of process_first_

  ASSERT("GroupProcessMpi::GroupProcessMpi",
	 "process_first out of range",
	 0 <= process_first_ && process_first_ < size);

  // Check range of process_last_plus_

  ASSERT("GroupProcessMpi::GroupProcessMpi",
	 "process_last_plus out of range",
	 process_last_plus_ == -1 || 
	 ( process_first_ < process_last_plus_ && 
	   process_last_plus_ <= size ));

  // Set default process_last_plus_ to P

  if (process_last_plus_ == -1) process_last_plus_ = size;

  // Compute size_

  size_ = process_last_plus_ - process_first_ ;

  int rank = Mpi::rank();

  rank_ = rank - process_first_;

}

//----------------------------------------------------------------------

void GroupProcessMpi::barrier()  const throw()
{ 
  MPI_Barrier (comm_); 
};

//----------------------------------------------------------------------

void GroupProcessMpi::sync (int rank, int tag) const throw()
{
  char buffer = 1;
  if (rank_ < rank) {
    call_mpi_ (__FILE__,__LINE__,"MPI_Send", MPI_Send 
		    (&buffer, 1, MPI_BYTE, rank, tag, comm_));
  } else if (rank_ > rank) {
    call_mpi_ (__FILE__,__LINE__,"MPI_Recv", MPI_Recv
		    (&buffer, 1, MPI_BYTE, rank, tag, comm_, MPI_STATUS_IGNORE));
  }
}

//----------------------------------------------------------------------

void * GroupProcessMpi::send_begin 
(int rank_dest, void * buffer, int size, int tag) const throw()
{
  MPI_Request * handle = 0;
  if (send_blocking_) {
    call_mpi_ (__FILE__,__LINE__,"MPI_Send",MPI_Send 
		    (buffer, size, MPI_BYTE, rank_dest, tag, comm_));
  } else {
    handle = new MPI_Request;
    call_mpi_ (__FILE__,__LINE__,"MPI_Isend",MPI_Isend
		    (buffer, size, MPI_BYTE, rank_dest, tag, comm_,handle));
  }
  return (void *)handle;
}

//----------------------------------------------------------------------

void GroupProcessMpi::send_end (void * handle) const throw()
{
  delete (MPI_Request *) handle;
}

//----------------------------------------------------------------------

void * GroupProcessMpi::recv_begin 
(int rank_source, void * buffer, int size, int tag) const throw()
{
  MPI_Request * handle = 0;
  if (recv_blocking_) {
    call_mpi_(__FILE__,__LINE__,"MPI_Recv",MPI_Recv
		   (buffer, size, MPI_BYTE, rank_source, tag, comm_, 
		    MPI_STATUS_IGNORE));
  } else {
    handle = new MPI_Request;
    call_mpi_(__FILE__,__LINE__,"MPI_Irecv", MPI_Irecv 
		   (buffer, size, MPI_BYTE, rank_source, tag, comm_, handle));
  }
  return handle;
}

//----------------------------------------------------------------------

void GroupProcessMpi::recv_end (void * handle) const throw()
{
  delete (MPI_Request *) handle;
}

//----------------------------------------------------------------------

void GroupProcessMpi::send_recv 
(int rank, void * buffer, int size, int tag) const throw()
{
  if (rank != process_rank_) {
    call_mpi_(__FILE__,__LINE__,"MPI_Sendrecv_replace",MPI_Sendrecv_replace
	      (buffer, size, MPI_BYTE, rank, tag, rank, tag,
	       comm_, MPI_STATUS_IGNORE));
  }
}

//----------------------------------------------------------------------

bool GroupProcessMpi::test (void * handle) const throw()
{
  int result = 1;
  if (handle != NULL) {
    call_mpi_(__FILE__,__LINE__,"MPI_Test",MPI_Test
		   ((MPI_Request*) (handle), &result, MPI_STATUS_IGNORE));
  }
  return result;
}

//----------------------------------------------------------------------

void GroupProcessMpi::wait (void * handle) const throw()
{
  if (handle != NULL) {
    call_mpi_(__FILE__,__LINE__,"MPI_Wait",MPI_Wait
		   ((MPI_Request*) (handle), MPI_STATUS_IGNORE));
  }
}

//----------------------------------------------------------------------

Reduce * GroupProcessMpi::create_reduce () const throw ()
{
  return new ReduceMpi (this);
}

//----------------------------------------------------------------------

void GroupProcessMpi::bulk_send_add
(int rank_dest, void * buffer, int size, int tag) const throw()
{
  INCOMPLETE("GroupProcessMpi::bulk_send_add");
}

//----------------------------------------------------------------------

void * GroupProcessMpi::bulk_send() const throw()
{
  INCOMPLETE("GroupProcessMpi::bulk_send");
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessMpi::bulk_recv_add
(int rank_source, void * buffer, int size, int tag) const throw()
{
  INCOMPLETE("GroupProcessMpi::bulk_recv_add");
}

//----------------------------------------------------------------------

void * GroupProcessMpi::bulk_recv() const throw()
{
  INCOMPLETE("GroupProcessMpi::bulk_recv");
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessMpi::bulk_wait(void * handle) const throw()
{
  INCOMPLETE("GroupProcessMpi::bulk_wait");
}

//======================================================================

void GroupProcessMpi::call_mpi_
(
 const char * file, 
 int          line , 
 const char * name, 
 int          ierr
 ) const throw()
  {
    if (ierr != MPI_SUCCESS) {
      char message[ERROR_LENGTH];
      char function[ERROR_LENGTH];
      sprintf (message,"MPI rank=%d  error=%d",rank_,ierr);
      sprintf (function,"GroupProcessMpi::%s",name);
      ERROR(function,message);
    }
  };

#endif /* CONFIG_USE_MPI */
