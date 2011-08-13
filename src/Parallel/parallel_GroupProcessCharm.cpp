// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcessCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 24 14:55:35 PDT 2011
/// @brief    Implementation of the GroupProcessCharm class


#ifdef CONFIG_USE_CHARM

//----------------------------------------------------------------------

#include "cello.hpp"

#include "parallel.hpp"

//----------------------------------------------------------------------

GroupProcessCharm::GroupProcessCharm
(int process_first,
 int process_last_plus) throw ()
  : GroupProcess(),
    process_first_    (process_first),
    process_last_plus_(process_last_plus)
{
  size_ = CkNumPes();
  rank_ = CkMyPe();
}

//----------------------------------------------------------------------

void GroupProcessCharm::barrier()  throw()
{ 
};

//----------------------------------------------------------------------

void GroupProcessCharm::sync (int rank, int tag) throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::send_begin 
(int rank_dest, void * buffer, int size, int tag) throw()
{
  return NULL;
}

//----------------------------------------------------------------------

void GroupProcessCharm::send_end (void * handle) throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::recv_begin 
(int rank_source, void * buffer, int size, int tag) throw()
{
  return NULL;
}

//----------------------------------------------------------------------

void GroupProcessCharm::recv_end (void * handle) throw()
{
}

//----------------------------------------------------------------------

bool GroupProcessCharm::test (void * handle) throw()
{
  return false;
}

//----------------------------------------------------------------------

void GroupProcessCharm::wait (void * handle) throw()
{
}

//----------------------------------------------------------------------

Reduce * GroupProcessCharm::create_reduce () throw ()
{
  return new ReduceCharm(this);
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_send_add(int rank_dest, void * buffer, int size, int tag) throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::bulk_send() throw()
{
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_recv_add(int rank_source, void * buffer, int size, int tag) throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::bulk_recv() throw()
{
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_wait(void * handle) throw()
{
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
