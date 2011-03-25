// $Id: parallel_GroupProcessCharm.cpp 2093 2011-03-12 01:17:05Z bordner $
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
  INCOMPLETE("GroupProcessCharm::GroupProcessCharm");
  size_ = CkNumPes();
  rank_ = CkMyPe();

}

//----------------------------------------------------------------------

void GroupProcessCharm::barrier()  throw()
{ 
  INCOMPLETE("GroupProcessCharm::barrier");
};

//----------------------------------------------------------------------

void GroupProcessCharm::sync (int rank, int tag) throw()
{
  INCOMPLETE("GroupProcessCharm::sync");
}

//----------------------------------------------------------------------

void * GroupProcessCharm::send_begin 
(int rank_dest, void * buffer, int size, int tag) throw()
{
  INCOMPLETE("GroupProcessCharm::send_begin");
  return NULL;
}

//----------------------------------------------------------------------

void GroupProcessCharm::send_end (void * handle) throw()
{
  INCOMPLETE("GroupProcessCharm::send_end");
}

//----------------------------------------------------------------------

void * GroupProcessCharm::recv_begin 
(int rank_source, void * buffer, int size, int tag) throw()
{
  INCOMPLETE("GroupProcessCharm::recv_begin");
  return NULL;
}

//----------------------------------------------------------------------

void GroupProcessCharm::recv_end (void * handle) throw()
{
  INCOMPLETE("GroupProcessCharm::recv_end");
}

//----------------------------------------------------------------------

bool GroupProcessCharm::test (void * handle) throw()
{
  INCOMPLETE("GroupProcessCharm::test");
  return false;
}

//----------------------------------------------------------------------

void GroupProcessCharm::wait (void * handle) throw()
{
  INCOMPLETE("GroupProcessCharm::wait");
}

//----------------------------------------------------------------------

Reduce * GroupProcessCharm::create_reduce () throw ()
{
  return new ReduceCharm(this);
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_send_add(int rank_dest, void * buffer, int size, int tag) throw()
{
  INCOMPLETE("GroupProcessCharm::bulk_send_add");
}

//----------------------------------------------------------------------

void * GroupProcessCharm::bulk_send() throw()
{
  INCOMPLETE("GroupProcessCharm::bulk_send");
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_recv_add(int rank_source, void * buffer, int size, int tag) throw()
{
  INCOMPLETE("GroupProcessCharm::bulk_recv_add");
}

//----------------------------------------------------------------------

void * GroupProcessCharm::bulk_recv() throw()
{
  INCOMPLETE("GroupProcessCharm::bulk_recv");
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_wait(void * handle) throw()
{
  INCOMPLETE("GroupProcessCharm::bulk_wait");
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
