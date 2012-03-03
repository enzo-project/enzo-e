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

void GroupProcessCharm::barrier()  const throw()
{ 
};

//----------------------------------------------------------------------

void GroupProcessCharm::sync (int rank, int tag) const throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::send_begin 
(int rank_dest, void * buffer, int size, int tag) const throw()
{
  return NULL;
}

//----------------------------------------------------------------------

void GroupProcessCharm::send_end (void * handle) const throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::recv_begin 
(int rank_source, void * buffer, int size, int tag) const throw()
{
  return NULL;
}

//----------------------------------------------------------------------

void GroupProcessCharm::recv_end (void * handle) const throw()
{
}

//----------------------------------------------------------------------

void GroupProcessCharm::send_recv (int rank, void * buffer, int size, int tag)
  const throw()
{
}

//----------------------------------------------------------------------

bool GroupProcessCharm::test (void * handle) const throw()
{
  return false;
}

//----------------------------------------------------------------------

void GroupProcessCharm::wait (void * handle) const throw()
{
}

//----------------------------------------------------------------------

Reduce * GroupProcessCharm::create_reduce () const throw ()
{
  return new ReduceCharm(this);
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_send_add
(int rank_dest, void * buffer, int size, int tag) const throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::bulk_send() const throw()
{
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_recv_add
(int rank_source, void * buffer, int size, int tag) const throw()
{
}

//----------------------------------------------------------------------

void * GroupProcessCharm::bulk_recv() const throw()
{
  return 0;
}

//----------------------------------------------------------------------

void GroupProcessCharm::bulk_wait(void * handle) const throw()
{
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
