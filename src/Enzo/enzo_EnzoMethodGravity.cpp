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

EnzoMethodGravity::EnzoMethodGravity
(const FieldDescr * field_descr,
 Solver * solver,
 double grav_const)
  : Method(),
    grav_const_(grav_const),
    solver_(NULL)
{
  const int num_fields = field_descr->field_count();
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields(num_fields);
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

  // Initialize the linear system

  Field field = block->data()->field();

  Matrix * A = new EnzoMatrixLaplace;
  const int ix = field.field_id ("potential");

  /// access problem-defining fields for eventual RHS and solution
  const int id  = field.field_id("density");
  const int idt = field.field_id("density_total");

  const int ib = field.field_id("potential");
  
    
  // Solve the linear system
  solver_->apply (A, ix, ib, block);
  
  block->compute_done();
}

