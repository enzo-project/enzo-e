// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Layout class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "parallel.hpp"

//----------------------------------------------------------------------

Layout::Layout() throw()
  : process_first_ (0),
    process_count_ (1)
{

  for (int i=0; i<3; i++) {
    block_count_[i] = 1;
  }
#ifdef CONFIG_USE_MPI
  mpi_comm_ = MPI_COMM_CELLO;
  MPI_Comm_group (mpi_comm_, &mpi_group_);
#endif
}

//----------------------------------------------------------------------

void Layout::set_process_range(int process_first, int process_count) throw()
{
  process_first_ = process_first;
  process_count_  = process_count;
}

//----------------------------------------------------------------------

void Layout::set_block_count(int nbx, int nby, int nbz) throw()
{
  block_count_[0] = nbx;
  block_count_[1] = nby;
  block_count_[2] = nbz; 
}

//----------------------------------------------------------------------

int Layout::block_count (int *nbx, int *nby, int *nbz) throw()
{ 
  *nbx = block_count_[0];
  *nby = block_count_[1];
  *nbz = block_count_[2];
  return (*nbx)*(*nby)*(*nbz);
}

//----------------------------------------------------------------------

void Layout::process_range(int * process_first, int * process_count) throw()
{
  *process_first = process_first_;
  *process_count  = process_count_;
}

//----------------------------------------------------------------------

int Layout::local_count (int ip) throw()
{
  int block_count = block_count_[0] * block_count_[1] * block_count_[2];
  if (process_first_ <= ip && ip < process_first_ + process_count_) {
    return (ip*block_count)/process_count_ - ((ip-1)*block_count)/process_count_;
  } else {
    return 0;
  }
  
}

//----------------------------------------------------------------------

int Layout::process (int ib)  throw()
{
  int block_count = block_count_[0] * block_count_[1] * block_count_[2];
  if (0 <= ib && ib < block_count) {
    return process_first_ + process_count_*ib / block_count;
  } else {
    return PROCESS_NULL;
  }
}

//----------------------------------------------------------------------

int Layout::block_index (int ibx, int iby, int ibz) throw()
{
  return ibx + block_count_[0]*(iby + block_count_[1]*ibz);
}

//----------------------------------------------------------------------

void Layout::block_indices (int ib, int * ibx, int * iby, int * ibz) throw()
{
  *ibx = ib % block_count_[0];
  *iby = (ib / block_count_[0]) % block_count_[1];
  *ibz = ib / (block_count_[0]*block_count_[1]);
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_MPI

void Layout::initialize_mpi_()
{

  UNTESTED_MESSAGE("Patch::set_layout");

  // Delete old communicator if needed

  if (mpi_comm_ != MPI_COMM_CELLO) {
    MPI_Comm_free  (&mpi_comm_);
    MPI_Group_free (&mpi_group_);
  }

  // Check if layout process range makes sense

  int layout_first, layout_size, mpi_size;

  process_range(&layout_first,&layout_size);
  MPI_Comm_size(MPI_COMM_CELLO, &mpi_size);

  printf ("%s:%d DEBUG %d %d %d\n",
	  __FILE__,__LINE__,layout_first, layout_size, mpi_size);
  
  if ( ! ((0 <= layout_first) &&
	  (layout_first + layout_size <= mpi_size))) {
    char buffer[ERROR_MESSAGE_LENGTH];
    sprintf (buffer,
	     "Illegal layout_first = %d layout_size = %d mpi_size = %d",
	     layout_first, layout_size, mpi_size);
  
    ERROR_MESSAGE("Patch::set_layout",  buffer);
  }

  if (layout_first == 0 && layout_size == mpi_size) {

    // Use MPI_COMM_CELLO for group / comm
    mpi_comm_ = MPI_COMM_CELLO;
    MPI_Comm_group (MPI_COMM_CELLO, &mpi_group_);

  } else {

    // Create new group / comm with layout's range of processes


    MPI_Group mpi_group;
    MPI_Comm_group (MPI_COMM_CELLO, &mpi_group);

    int ranges[1][3] = {{layout_first,layout_size - layout_first,1}}
;
    MPI_Group_range_incl(mpi_group,1,ranges,&mpi_group_);
    MPI_Comm_create (MPI_COMM_CELLO, mpi_group_, &mpi_comm_);
  }
}
#endif
