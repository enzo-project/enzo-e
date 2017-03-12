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

// #define DEBUG_METHOD

#ifdef DEBUG_METHOD
#   define TRACE_METHOD(METHOD,BLOCK)					\
  CkPrintf ("%d %s:%d %s TRACE %s %p\n",CkMyPe(),__FILE__,__LINE__, \
	    BLOCK->name().c_str(),METHOD,this);			    \
  fflush(stdout);
#else
#   define TRACE_METHOD(METHOD,BLOCK) /*  */ 
#endif
//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity
(const FieldDescr * field_descr,
 int index_solver,
 double grav_const)
  : Method(),
    index_solver_(index_solver),
    grav_const_(grav_const)
{
  const int num_fields = field_descr->field_count();
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields(num_fields);
}

//----------------------------------------------------------------------

EnzoMethodGravity::~EnzoMethodGravity() throw()
{
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{

  TRACE_METHOD("compute()",block);

  // Initialize the linear system

  Field field = block->data()->field();

  Matrix * A = new EnzoMatrixLaplace;
  const int ix = field.field_id ("potential");

  /// access problem-defining fields for eventual RHS and solution
  const int id  = field.field_id("density");
  const int idt = field.field_id("density_total");
  const int idensity = (idt != -1) ? idt : id;
  const int ib = field.field_id ("B");
  
  // Solve the linear system
  field.scale(ib, -4.0 * (cello::pi) * grav_const_, idensity);

  Solver * solver = block->simulation()->problem()->solver(index_solver_);
  
  // May exit before solve is done...
  solver->set_callback (CkIndex_EnzoBlock::r_method_gravity_continue(NULL));

  solver->apply (A, ix, ib, block);

}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_continue(CkReductionMsg * msg)
{

  TRACE_METHOD("r_method_gravity_end()",this);
  
  // So do refresh with barrier synch (note barrier instead of
  // neighbor synchronization otherwise will conflict with Method
  // refresh ("Charm++ fatal error: mis-matched client callbacks in
  // reduction messages")

  Refresh refresh (4,0,neighbor_leaf, sync_barrier);
  refresh.set_active(is_leaf());
  refresh.add_all_fields(data()->field().field_count());

  refresh_enter(CkIndex_EnzoBlock::r_method_gravity_end(NULL),&refresh);

  //  delete msg;
}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_end(CkReductionMsg * msg)
{
  TRACE_METHOD("r_method_gravity_end()",this);
  
  delete msg;
  
  // BUG: acceleration computed before Solver completes
  
  /// compute acceleration fields from potential
  int order;
  EnzoComputeAcceleration compute_acceleration(data()->field().field_descr(),
					       rank(), order=4);

  compute_acceleration.compute(this);

  // wait for all Blocks before continuing
  compute_done();
}

