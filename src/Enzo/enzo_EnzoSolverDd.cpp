// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDd.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-10-01
/// @brief    Implements the EnzoSolverDd class
///
/// @brief [\ref Enzo] Multigrid on the root-level grid using Dd, then
/// BiCgStab in overlapping subdomains defined by root-level Blocks.
/// An optional final Jacobi step can be applied to smooth the solution
/// along subdomain boundaries.
///
/// 0. restrict b to level 0
/// 1. coarse solve: solve Ax=b on level 0
/// 2. Prolong x to child blocks as xc
/// 3. domain solve: solve Ai xi = bi in each root-grid block
///        use xc for boundary conditions and initial guess
/// 4. final smoother: apply final smoother on A x = b (if any)
///
///

#include "cello.hpp"
#include "enzo.hpp"

//======================================================================

EnzoSolverDd::EnzoSolverDd
  (std::string name,
   std::string field_x,
   std::string field_b,
   int monitor_iter,
   int restart_cycle,
   int solve_type,
   int index_prolong,
   int index_restrict,
   int min_level,
   int max_level,
   int index_solve_coarse,
   int index_solve_domain,
   int index_solve_smooth,
   int coarse_level)
    : Solver(name,
	     field_x,
	     field_b,
	     monitor_iter,
	     restart_cycle,
	     solve_type,
             index_prolong,
             index_restrict,
	     min_level,
	     max_level),
      index_solve_coarse_(index_solve_coarse),
      index_solve_domain_(index_solve_domain),
      index_solve_smooth_(index_solve_smooth),
      index_prolong_(index_prolong),
      index_restrict_(index_restrict),
      ixc_(-1),
      mx_(0),my_(0),mz_(0),
      gx_(0),gy_(0),gz_(0),
      coarse_level_(coarse_level)
{

  ixc_ = cello::define_field("X_c");

  /// Initialize default Refresh

  Refresh * refresh = cello::refresh(ir_post_);
  cello::simulation()->refresh_set_name(ir_post_,name);

  refresh->add_field (ix_);

  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  i_sync_restrict_ = scalar_descr_sync->new_value(name + ":restrict");
  i_sync_prolong_  = scalar_descr_sync->new_value(name + ":prolong");

  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  i_msg_prolong_ = scalar_descr_void->new_value(name + ":msg_prolong");
  for (int ic=0; ic<cello::num_children(); ic++) {
    i_msg_restrict_[ic] = scalar_descr_void->new_value(name + ":msg_restrict");
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::apply ( std::shared_ptr<Matrix> A, Block * block) throw()
{
  Solver::begin_(block);

  Field field = block->data()->field();

  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  A_ = A;

  allocate_temporary_(block);

  // Check that component solvers are of the correct type
  ASSERT2("EnzoSolverDd::apply()",
	  "Coarse solver %s type %s != solve_level",
	  cello::solver(index_solve_coarse_)->name().c_str(),
	  solve_string[cello::solver(index_solve_coarse_)->solve_type()],
	  cello::solver(index_solve_coarse_)->solve_type() == solve_level);

  Sync * sync_restrict = psync_restrict(block);

  sync_restrict->set_stop(1 + cello::num_children()); // self and children

  Sync * sync_prolong = psync_prolong(block);

  sync_prolong->set_stop(1 + 1); // self and parent

  int level = block->level();

  const int m = mx_*my_*mz_;

  if ( ! block->is_leaf() ) {
    std::fill_n ((enzo_float*) field.values(ib_), m, 0.0);
  }

  std::fill_n ((enzo_float*) field.values(ix_),  m, 0.0);
  std::fill_n ((enzo_float*) field.values(ixc_), m, 0.0);
	
  if (block->is_leaf()) {

    begin_solve(enzo::block(block));

  } else {

    const bool in_range = (coarse_level_ <= level && level <= max_level_);

    if ( ! in_range ) {

      call_coarse_solver(enzo::block(block));

    } else {

      // non-leaf blocks in range [coarse_level,max_level_]

      // count self for restrict receive
      restrict_recv(enzo::block(block),nullptr);
    }
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::begin_solve(EnzoBlock * enzo_block) throw()
{
  if (enzo_block->level() == coarse_level_) {

    call_coarse_solver(enzo_block);

  } else {

    restrict (enzo_block);

  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::restrict(EnzoBlock * enzo_block) throw()
{
  restrict_send(enzo_block);

  call_coarse_solver(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverDd::restrict_send(EnzoBlock * enzo_block) throw()
{
  // Pack field
  Index index = enzo_block->index();
  int level   = index.level();
  int ic3[3];
  index.child(level,&ic3[0],&ic3[1],&ic3[2],min_level_);

  FieldMsg * msg = pack_field_(enzo_block,ib_,refresh_coarse,ic3);

  // Send packed field to parent
  Index index_parent = enzo_block->index().index_parent(min_level_);
  enzo::block_array()[index_parent].p_solver_dd_restrict_recv(msg);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_restrict_recv(FieldMsg * msg)
{
  static_cast<EnzoSolverDd*> (solver())->restrict_recv(this,msg);
}

//----------------------------------------------------------------------

void EnzoSolverDd::restrict_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
{
  // Unpack "B" vector data from children

  // Save field message from child
  if (msg != NULL) *pmsg_restrict(enzo_block,msg->child_index()) = msg;

  // Continue if all expected messages received
  if (psync_restrict(enzo_block)->next() ) {

    // Restore saved messages then clear
    for (int i=0; i<cello::num_children(); i++) {
      msg = *pmsg_restrict(enzo_block,i);
      *pmsg_restrict(enzo_block,i) = NULL;
      // Unpack field from message then delete message
      unpack_field_(enzo_block,msg,ib_,refresh_coarse);
    }

    begin_solve(enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::call_coarse_solver(EnzoBlock * enzo_block) throw()
{
  Solver * solve_coarse = cello::solver(index_solve_coarse_);

  solve_coarse->set_sync_id (enzo_sync_id_solver_dd_coarse);
  solve_coarse->set_callback(CkIndex_EnzoBlock::p_solver_dd_solve_coarse());

  solve_coarse->set_field_x (ixc_);
  solve_coarse->set_field_b (ib_);

  solve_coarse->apply(A_,enzo_block);
}


//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_solve_coarse()
{
  CkCallback callback(CkIndex_EnzoBlock::r_solver_dd_barrier(NULL),
		      enzo::block_array());
  contribute(callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_dd_barrier(CkReductionMsg * msg)
{
  static_cast<EnzoSolverDd*> (solver())->prolong(this);
  delete msg;
}

//----------------------------------------------------------------------

void EnzoSolverDd::prolong(EnzoBlock * enzo_block) throw()
{

  if (is_finest_(enzo_block)) {
    copy_xc_to_x_(enzo_block);
  }

  /// Prolong solution to next-finer level

  const int level = enzo_block->level();

  if (level == coarse_level_) {

    prolong_send_ (enzo_block);

  }

  if (coarse_level_ < level && level <= max_level_) {
    enzo_block->solver_dd_prolong_recv(NULL);
  } else {
    call_domain_solver (enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::prolong_send_(EnzoBlock * enzo_block) throw()
{
  if ( ! is_finest_(enzo_block) ) {
    ItChild it_child(cello::rank());
    int ic3[3];

    while (it_child.next(ic3)) {

      FieldMsg * msg = pack_field_(enzo_block,ixc_,refresh_fine,ic3);

      Index index_child = enzo_block->index().index_child(ic3,min_level_);

      enzo::block_array()[index_child].p_solver_dd_prolong_recv(msg);

    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_prolong_recv(FieldMsg * msg)
{  solver_dd_prolong_recv(msg); }

void EnzoBlock::solver_dd_prolong_recv(FieldMsg * msg)
{
  static_cast<EnzoSolverDd*> (solver())->prolong_recv(this,msg);
}

//----------------------------------------------------------------------

void EnzoSolverDd::prolong_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
{
  // Unpack "X" vector data from parent

  // Save field message from parent
  if (msg != NULL) *pmsg_prolong(enzo_block) = msg;

  // Continue if all expected messages received
  if (psync_prolong(enzo_block)->next() ) {

    // Restore saved message then clear
    msg = *pmsg_prolong(enzo_block);
    *pmsg_prolong(enzo_block) = NULL;

    // Unpack field from message then delete message
    unpack_field_(enzo_block,msg,ixc_,refresh_fine);

    // copy X = XC
    // copy X_copy = XC (using Solver::reuse_solution_(cycle) )
    copy_xc_to_x_(enzo_block);

    call_domain_solver(enzo_block);

    prolong_send_ (enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::copy_xc_to_x_(EnzoBlock * enzo_block) throw()
{

  const int m = mx_*my_*mz_;

  Field field = enzo_block->data()->field();

  std::copy_n((enzo_float *) field.values(ixc_),m,
	      (enzo_float *) field.values(ix_));
  std::copy_n((enzo_float *) field.values(ixc_),m,
	      (enzo_float *) field.values("X_copy"));

}

//----------------------------------------------------------------------

void EnzoSolverDd::call_domain_solver(EnzoBlock * enzo_block) throw()
{
  Solver * solve_domain = cello::solver(index_solve_domain_);

  solve_domain->set_sync_id (enzo_sync_id_solver_dd_domain);
  solve_domain->set_callback(CkIndex_EnzoBlock::p_solver_dd_solve_domain());

  solve_domain->set_field_x (ix_);
  solve_domain->set_field_b (ib_);

  solve_domain->apply(A_,enzo_block);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_solve_domain()
{
  static_cast<EnzoSolverDd*> (solver())->continue_after_domain_solve(this);
}

//----------------------------------------------------------------------

void EnzoSolverDd::continue_after_domain_solve(EnzoBlock * enzo_block) throw()
{
  CkCallback callback(CkIndex_EnzoBlock::r_solver_dd_end(NULL),
		      enzo::block_array());
  enzo_block->contribute(callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_dd_end(CkReductionMsg * msg)
{
  static_cast<EnzoSolverDd*> (solver())->call_last_smoother(this);
  delete msg;
}

//----------------------------------------------------------------------

void EnzoSolverDd::call_last_smoother(EnzoBlock * enzo_block) throw()
{
  Solver * smooth_last = cello::solver(index_solve_smooth_);

  smooth_last->set_sync_id (enzo_sync_id_solver_dd_smooth);
  smooth_last->set_callback(CkIndex_EnzoBlock::p_solver_dd_last_smooth());

  smooth_last->set_field_x(ix_);
  smooth_last->set_field_b(ib_);

  smooth_last->apply(A_,enzo_block);

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_last_smooth()
{
  static_cast<EnzoSolverDd*> (solver())->continue_after_last_smooth(this);
}

//----------------------------------------------------------------------

void EnzoSolverDd::continue_after_last_smooth(EnzoBlock * enzo_block) throw()
{
  end(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverDd::end (Block* block) throw ()
{
  deallocate_temporary_(block);

  Solver::end_(block);
}

//======================================================================

FieldMsg * EnzoSolverDd::pack_field_(EnzoBlock * enzo_block,
				     int index_field,
				     int refresh_type,
				     int * ic3)
{
  int  if3[3] = {0,0,0};
  int g3[3];
  cello::field_descr()->ghost_depth(index_field,g3,g3+1,g3+2);
  if (refresh_type != refresh_fine)
    for (int i=0; i<3; i++) g3[i]=0;

  Refresh * refresh = new Refresh;
  refresh->set_prolong(index_prolong_);
  refresh->set_restrict(index_restrict_);
  refresh->add_field(index_field);

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, g3, refresh_type, refresh, true);

  if (refresh_type == refresh_fine) {
    refresh->set_prolong(index_prolong_);
  } else if (refresh_type == refresh_coarse) {
    refresh->set_restrict(index_restrict_);
  }

  Field field = enzo_block->data()->field();
  int narray;
  char * array;
  field_face->face_to_array(field,&narray,&array);

  delete field_face;

  FieldMsg * msg  = new (narray) FieldMsg;

  msg->n = narray;
  memcpy (msg->a, array, narray);
  delete [] array;

  msg->ic3[0] = ic3[0];
  msg->ic3[1] = ic3[1];
  msg->ic3[2] = ic3[2];

  return msg;

}

//----------------------------------------------------------------------

void EnzoSolverDd::unpack_field_
(EnzoBlock * enzo_block,
 FieldMsg * msg,
 int index_field,
 int refresh_type)
{
  int if3[3] = {0,0,0};
  int g3[3];
  cello::field_descr()->ghost_depth(index_field,g3,g3+1,g3+2);
  if (refresh_type != refresh_fine)
    for (int i=0; i<3; i++) g3[i]=0;
  Refresh * refresh = new Refresh;
  refresh->set_prolong(index_prolong_);
  refresh->set_restrict(index_restrict_);
  refresh->add_field(index_field);

  int * ic3 = msg->ic3;

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, g3, refresh_type, refresh, true);

  if (refresh_type == refresh_fine) {
    refresh->set_prolong(index_prolong_);
  } else if (refresh_type == refresh_coarse) {
    refresh->set_restrict(index_restrict_);
  }

  Field field = enzo_block->data()->field();

  char * a = msg->a;
  field_face->array_to_face(a, field);
  delete field_face;

  delete msg;
}

