// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodWriteCheck.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2202-02-12
/// @brief    Implements the EnzoMethodWriteCheck class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodWriteCheck::EnzoMethodWriteCheck ()
  : Method()
{
}

//----------------------------------------------------------------------

void EnzoMethodWriteCheck::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodWriteCheck::compute ( Block * block) throw()
{
  CkPrintf ("TRACE_METHOD_WRITE_CHECK %s compute()\n",block->name().c_str());

  block->compute_done();

}

