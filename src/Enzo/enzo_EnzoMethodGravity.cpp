// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class


#include "cello.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#ifdef DEBUG_METHOD
#   define TRACE_METHOD(method)						\
  CkPrintf ("%d %s:%d TRACE %s %p\n",CkMyPe(),__FILE__,__LINE__,method,this); \
  fflush(stdout);
#else
#   define TRACE_METHOD(method) /*  */ 
#endif
//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity(Solver * solver)
  : Method(),
    solver_(NULL)
{
}

//----------------------------------------------------------------------

EnzoMethodGravity::~EnzoMethodGravity() throw()
{
  delete solver_;
  solver_ = NULL;
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{
  block->compute_done();
}

