// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcessSerial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Mar  4 17:15:36 PST 2011
/// @brief    Implementation of the GroupProcessSerial class

#include "cello.hpp"

#include "parallel.hpp"

//----------------------------------------------------------------------

void * GroupProcessSerial::send_begin
(
 int    rank_dest, 
 void * buffer, 
 int    size, 
 int    tag
 ) throw()
{
  if (buffer_[tag] != 0) {
    WARNING("send_begin",
	    "multiple sends with no corresponding receive");
  }
  buffer_[(long int)tag] = buffer;
  return (void * ) tag;
}

//----------------------------------------------------------------------

void GroupProcessSerial::recv_end(void * handle) throw()
{
  if (buffer_[(long int)(handle)] == 0) {
    WARNING("recv_end",
	    "receive with no corresponding send");
  }
}

//----------------------------------------------------------------------

void GroupProcessSerial::send_recv (int rank, void * buffer, int size, int tag)
  throw()
{
  return;
}

//----------------------------------------------------------------------

Reduce * GroupProcessSerial::create_reduce () throw ()
{
  return new ReduceSerial (this);
}
