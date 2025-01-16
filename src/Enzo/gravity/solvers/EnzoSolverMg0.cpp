// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverMg0.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoSolverMg0 class
///
/// Multigrid solver on a non-adaptive mesh.  Can be any mesh level, but
/// typically the root-grid (level = 0).
///
///======================================================================
///
///  "Coarse" view of MG0 multigrid solver
///
///   @code
///
///   $MG(A_h,X_h,B_h)$
///
///    while ( ! converged() )
///       if (level == min_level) then
///          solve_coarse()     solve $A_h X_h = B_h$
///       else
/// 1        p_pre_smooth()     smooth $A_h X_h = B_h$
/// 2        p_residual()       $R_h = B_h - A_h * X_h$
/// 3        p_restrict ()      $B_H = I_h^H R_h$
/// 4        MG()               solve $A_H X_H = B_H$  (repeat for W-cycle)
/// 5        p_prolong ()       $X_h = X_h + I_H^h X_H$
/// 6        p_post_smooth()    smooth $A_h X_h = B_h$
///
///  @endcode
///
///----------------------------------------------------------------------
///
///  "Fine" view of MG0 multigrid solver
///
///  @code
///
///  enter_solver()
///
///     iter = 0
///     initialize X,R,C
///     if (level == max_level)
///        begin_cycle()
///
///  begin_cycle()
///
///     if (converged()) exit()
///     if (level == min_level) then
///        solve_coarse(A,X,B)
///     else
///        callback = p_pre_smooth()
///        call refresh (X,"level")
///
///  p_pre_smooth()
///
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
///
///  p_restrict_send(X)
///
///      A.residual(R,B,X) on level
///      pack R
///      index_parent.p_restrict_recv(R)
///
///  p_restrict_recv(B)
///
///      unpack B
///      --level
///      if (sync_restrict.next())
///          begin_cycle()
///
///  coarse_solve(A,X,B)
///
///      solve A X = B
///      prolong_send(X)
///
///  prolong_send(X)
///
///      if (level < max_level)
///         for child
///            pack X
///            child.prolong_recv(X)
///      else
///         begin_cycle()
///
///  prolong_recv(C)
///
///      ++level
///      unpack C
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
///
///  p_post_smooth(A,X,B)
///
///      smooth.apply (A,X,B)
///      prolong_send()
///
///  @endcode
///
///======================================================================

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/gravity/gravity.hpp"

//======================================================================

EnzoSolverMg0::EnzoSolverMg0
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
 int iter_max,
 double res_tol,
 int index_smooth_pre,
 int index_solve_coarse,
 int index_smooth_post,
 int index_smooth_last,
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
    bs_(0), bc_(0),
    rr_(0), rr_local_(0), rr0_(0),
    res_tol_(res_tol),
    A_(nullptr),
    index_smooth_pre_(index_smooth_pre),
    index_solve_coarse_(index_solve_coarse),
    index_smooth_post_(index_smooth_post),
    index_smooth_last_(index_smooth_last),
    iter_max_(iter_max),
    ic_(-1), ir_(-1),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    coarse_level_(coarse_level)
{
  // Initialize temporary fields

  ir_ = cello::field_descr()->insert_temporary();
  ic_ = cello::field_descr()->insert_temporary();

  Refresh * refresh = cello::refresh(ir_post_);
  cello::simulation()->refresh_set_name(ir_post_,name);

  refresh->add_field (ix_);
  refresh->add_field (ir_);
  refresh->add_field (ic_);

  refresh->set_min_face_rank(cello::rank() - 1);
  
  ScalarDescr * scalar_descr_int  = cello::scalar_descr_int();
  i_iter_  = scalar_descr_int ->new_value(name + ":iter");

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

void EnzoSolverMg0::apply ( std::shared_ptr<Matrix> A, Block * block) throw()
{
  Solver::begin_(block);

  A_ = A;

  allocate_temporary_(block);

  // clear scalars
  bs_ = 0.0;
  bc_ = 0.0;
  rr_ = 0.0;
  rr_local_ = 0.0;
  rr0_ = 0.0;
  *piter(block) = 0.0;

  /// Current and initial residual norm R'*R

  Field field = block->data()->field();

  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  EnzoBlock * enzo_block = enzo::block(block);

  // Initialize sync counters for restrict and prolong

  Sync * sync_restrict = psync_restrict(block);

  sync_restrict->set_stop(1 + cello::num_children()); // self and children

  Sync * sync_prolong = psync_prolong(block);

  sync_prolong->set_stop(1 + 1); // self and parent
  enter_solver_ (enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::enter_solver_ (EnzoBlock * enzo_block) throw()
///     iter = 0
///     initialize X,B,R,C
///     if (level == max_level)
///        begin_cycle()
{
  *piter(enzo_block) = 0.0;

  Field field = enzo_block->data()->field();

  enzo_float * X = (enzo_float*) field.values(ix_);
  enzo_float * R = (enzo_float*) field.values(ir_);
  enzo_float * C = (enzo_float*) field.values(ic_);

  // X = 0
  // R = B ( residual with X = 0 )
  // C = 0

  std::fill_n(X,mx_*my_*mz_,0.0);
  std::fill_n(R,mx_*my_*mz_,0.0);
  std::fill_n(C,mx_*my_*mz_,0.0);

  if (A_->is_singular()) {

    // Compute sum(B) and length() to project B onto range of A
    // if A is singular (

    long double reduce[2] = {0.0, 0.0};

    if (is_finest_(enzo_block)) {

      compute_shift_(enzo_block,reduce);
    }

    /// initiate callback for p_solver_begin_solve and contribute to
    /// sum and count

    CkCallback callback(CkIndex_EnzoBlock::r_solver_mg0_begin_solve(nullptr),
                        enzo::block_array());

    PERF_REDUCE_START(perf_reduce_solver_mg0);
    enzo_block->contribute(2*sizeof(long double), &reduce,
                           sum_long_double_2_type, callback);
  } else {

    begin_solve (enzo_block,nullptr);

  }

}

//----------------------------------------------------------------------

void EnzoSolverMg0::compute_shift_
(EnzoBlock * enzo_block,long double * reduce) throw()
{
  Field field = enzo_block->data()->field();

  enzo_float* B = (enzo_float*) field.values(ib_);

  for (int iz=gz_; iz<mz_-gz_; iz++) {
    for (int iy=gy_; iy<my_-gy_; iy++) {
      for (int ix=gx_; ix<mx_-gx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	reduce[0] += B[i];
	reduce[1] += 1.0;
      }
    }
  }
}
//----------------------------------------------------------------------

void EnzoBlock::r_solver_mg0_begin_solve(CkReductionMsg* msg)
{
  PERF_REDUCE_STOP(perf_reduce_solver_mg0);
  static_cast<EnzoSolverMg0*> (solver())->begin_solve(this,msg);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::begin_solve(EnzoBlock * enzo_block,
				CkReductionMsg *msg) throw()
{
  do_shift_(enzo_block,msg);

  // control flow starts at leaves, even in level > max_level,
  // since coarse solve may require reductions over all Blocks

  if (is_finest_(enzo_block)) {

    begin_cycle_ (enzo_block);

  } else {


    int level = enzo_block->level();
    if ( coarse_level_ <= level && level <= max_level_) {
      restrict_recv(enzo_block,nullptr);
    } else {
      call_coarse_solver(enzo_block);
    }

  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::do_shift_(EnzoBlock * enzo_block,
			      CkReductionMsg *msg) throw()
{
  if (msg != nullptr) {

    long double* data = (long double*) msg->getData();

    bs_ = data[0];
    bc_ = data[1];

    delete msg;
  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::begin_cycle_(EnzoBlock * enzo_block) throw()
///     if (converged()) exit()
///     if (level == min_level) then
///        coarse_solve(A,X,B)
///     else
///        callback = p_pre_smooth()
///        call refresh (X,"level")
{
  monitor_output_(enzo_block);

  Field field = enzo_block->data()->field();

  const int level = enzo_block->level();

  if ( ! is_finest_(enzo_block) ) {
    enzo_float * X = (enzo_float*) field.values(ix_);
    std::fill_n(X,mx_*my_*mz_,0.0);
  }

  if (level == coarse_level_) {

    call_coarse_solver(enzo_block);

  } else {

    if (index_smooth_pre_ >= 0) {

      call_pre_smoother (enzo_block);

    } else {

      restrict (enzo_block);

    }

  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::monitor_output_(EnzoBlock * enzo_block)
{
  const int iter = *(piter(enzo_block));

  const bool l_output =
    ( ( enzo_block->index().is_root()) &&
      ( (iter == 0))); // ||

  if (l_output) {
    Solver::monitor_output_(enzo_block,iter,rr0_,0.0,rr_,0.0);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_solve_coarse()
{
  EnzoSolverMg0 * solver =
    static_cast<EnzoSolverMg0*> (this->solver());

  CkCallback callback(CkIndex_EnzoBlock::r_solver_mg0_barrier(nullptr),
		      enzo::block_array());

  long double data[1] = {solver->rr_local()};

  PERF_REDUCE_START(perf_reduce_solver_mg0);
  contribute(sizeof(long double), data,  sum_long_double_type, callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_mg0_barrier(CkReductionMsg* msg)
{
  PERF_REDUCE_STOP(perf_reduce_solver_mg0);
  EnzoSolverMg0 * solver =
    static_cast<EnzoSolverMg0*> (this->solver());

  long double rr = ((long double*) msg->getData())[0];
  solver->set_rr(rr);
  solver->set_rr_local(0.0);
  if (*solver->piter(this)==0) solver->set_rr0(rr);

  delete msg;

  solver->prolong(this);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_restrict()
{
  EnzoSolverMg0 * solver =
    static_cast<EnzoSolverMg0*> (this->solver());

  solver->restrict(this);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::restrict(EnzoBlock * enzo_block) throw()
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  restrict_send (enzo_block);

  // All Blocks must call coarse solver since may involve
  // global reductions

  call_coarse_solver(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_coarse_solver(EnzoBlock * enzo_block) throw()
{
  Solver * solve_coarse = cello::solver(index_solve_coarse_);

  solve_coarse->set_min_level(min_level_);
  solve_coarse->set_max_level(coarse_level_);
  solve_coarse->set_sync_id (enzo_sync_id_solver_mg0_coarse);
  solve_coarse->set_callback(CkIndex_EnzoBlock::p_solver_mg0_solve_coarse());

  solve_coarse->set_field_x (ix_);
  solve_coarse->set_field_b (ib_);

  solve_coarse->apply(A_,enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_pre_smoother(EnzoBlock * enzo_block) throw()
{
  Solver * smooth_pre = cello::solver(index_smooth_pre_);

  smooth_pre->set_min_level(enzo_block->level());
  smooth_pre->set_max_level(enzo_block->level());
  smooth_pre->set_sync_id (enzo_sync_id_solver_mg0_pre);
  smooth_pre->set_callback(CkIndex_EnzoBlock::p_solver_mg0_restrict());

  smooth_pre->set_field_x(ix_);
  smooth_pre->set_field_b(ib_);

  smooth_pre->apply(A_,enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_post_smoother(EnzoBlock * enzo_block) throw()
{
  Solver * smooth_post = cello::solver(index_smooth_post_);

  smooth_post->set_min_level(enzo_block->level());
  smooth_post->set_max_level(enzo_block->level());
  smooth_post->set_sync_id (enzo_sync_id_solver_mg0_post);
  smooth_post->set_callback(CkIndex_EnzoBlock::p_solver_mg0_post_smooth());

  smooth_post->set_field_x(ix_);
  smooth_post->set_field_b(ib_);

  smooth_post->apply(A_,enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_last_smoother(EnzoBlock * enzo_block) throw()
{
  Solver * smooth_last = cello::solver(index_smooth_last_);

  smooth_last->set_sync_id (enzo_sync_id_solver_mg0_last);
  smooth_last->set_callback(CkIndex_EnzoBlock::p_solver_mg0_last_smooth());

  smooth_last->set_field_x(ix_);
  smooth_last->set_field_b(ib_);

  smooth_last->apply(A_,enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::restrict_send(EnzoBlock * enzo_block) throw()
///
///      A.residual(R,B,X)
///      pack R
///      index_parent.p_restrict_recv(R)
{
  compute_residual_(enzo_block);

  FieldMsg * msg = pack_residual_(enzo_block);

  Index index_parent = enzo_block->index().index_parent(min_level_);

  enzo::block_array()[index_parent].p_solver_mg0_restrict_recv(msg);

}

//----------------------------------------------------------------------

void EnzoSolverMg0::compute_residual_(EnzoBlock * enzo_block) throw()
{
  Field field = enzo_block->data()->field();

  A_->residual(ir_, ib_, ix_, enzo_block);

  if ( is_finest_(enzo_block) ) {
    enzo_float * R = (enzo_float*) field.values(ir_);
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  rr_local_ += R[i]*R[i];
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_restrict_recv(FieldMsg * msg)
{
  EnzoSolverMg0 * solver =
    static_cast<EnzoSolverMg0*> (this->solver());

  solver->restrict_recv(this,msg);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::restrict_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
///
///      [ unpack B ]
///      if (sync.next())
///          begin_cycle()
{

  // Unpack "B" vector data from children

  // Save field message from child
  if (msg != nullptr) *pmsg_restrict(enzo_block,msg->child_index()) = msg;

  // Continue if all expected messages received
  if (psync_restrict(enzo_block)->next() ) {

    // Restore saved messages then clear
    for (int i=0; i<cello::num_children(); i++) {
      msg = *pmsg_restrict(enzo_block,i);
      *pmsg_restrict(enzo_block,i) = nullptr;
      // Unpack field from message then delete message
      unpack_residual_(enzo_block,msg);
    }

    begin_cycle_ (enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::prolong(EnzoBlock * enzo_block) throw()
///
///      solve A X = B
///      end_cycle()
{
  /// Prolong solution to next-finer level

  const int level = enzo_block->level();

  if (level == coarse_level_) {

    if ( ! is_finest_(enzo_block) ) {

      prolong_send_ (enzo_block);

    }
  }

  if (coarse_level_ < level && level <= max_level_) {
    prolong_recv(enzo_block,nullptr);
  } else {

    end_cycle (enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::prolong_send_(EnzoBlock * enzo_block) throw()
///
///      for child
///         pack X
///         child.prolong_recv(X)
{
  ItChild it_child(cello::rank());
  int ic3[3];

  while (it_child.next(ic3)) {

    FieldMsg * msg = pack_correction_(enzo_block,ic3);

    Index index_child = enzo_block->index().index_child(ic3,min_level_);

    enzo::block_array()[index_child].p_solver_mg0_prolong_recv(msg);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg * msg)
{
  static_cast<EnzoSolverMg0*> (solver())->prolong_recv(this,msg);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::prolong_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
///
///      [ unpack C ]
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{

  // Save message

  // Return if not ready yet
  if (msg != nullptr) *pmsg_prolong(enzo_block) = msg;

  if (! psync_prolong(enzo_block)->next() ) return;

  // Restore saved message then clear
  msg = *pmsg_prolong(enzo_block);
  *pmsg_prolong(enzo_block) = nullptr;
  // Unpack "C" vector data from children
  unpack_correction_(enzo_block,msg);

  Field field = enzo_block->data()->field();

  enzo_float * X = (enzo_float*) field.values(ix_);
  enzo_float * C = (enzo_float*) field.values(ic_);

  for (int i=0; i<mx_*my_*mz_; i++) {
    X[i] += C[i];
  }

  if (index_smooth_post_ >= 0) {

    call_post_smoother(enzo_block);

  } else {

    post_smooth (enzo_block);
  }

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_post_smooth()
{
  EnzoSolverMg0 * solver =
    static_cast<EnzoSolverMg0*> (this->solver());

  solver->post_smooth(this);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::post_smooth(EnzoBlock * enzo_block) throw()
///
///      smooth.apply (A,X,B)
///      end_cycle()
{
  if ( ! is_finest_(enzo_block) ) {

    prolong_send_ (enzo_block);
  }

  end_cycle (enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::end_cycle(EnzoBlock * enzo_block) throw()
///
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
{
  ++ (*piter(enzo_block));

  bool is_converged = is_converged_(enzo_block);
  bool is_diverged  = is_diverged_(enzo_block);

  const int iter = *piter(enzo_block);
	
  const bool l_output =
    ( ( enzo_block->index().is_root()) &&
      ( (is_converged) || (is_diverged) ||
	(monitor_iter_ && (iter % monitor_iter_) == 0 )) );

  if (l_output) {
    Solver::monitor_output_(enzo_block,iter,rr0_,0.0,rr_,0.0);
  }

  if (is_converged || is_diverged) {

    // Do an optional final smoothing on the full mesh For use in Dan
    // Reynolds HG algorithm in which Mg0 with no pre- or
    // post-smoothings is used as a preconditioner to BiCgStab

    if (index_smooth_last_ >= 0 && (is_finest_(enzo_block)) ) {

      call_last_smoother(enzo_block);

    } else {

      end (enzo_block);
    }

  } else {

    if ( is_finest_(enzo_block)) {

      begin_cycle_ (enzo_block);

    } else {

      int level = enzo_block->level();
      if ( ! (coarse_level_ <= level && level <= max_level_) ) {
	call_coarse_solver(enzo_block);
      } else {
        restrict_recv(enzo_block,nullptr);
      }

    }

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_last_smooth()
{
  EnzoSolverMg0 * solver = static_cast<EnzoSolverMg0*> (this->solver());

  solver->end(this);
}

//----------------------------------------------------------------------

FieldMsg * EnzoSolverMg0::pack_residual_(EnzoBlock * enzo_block) throw()
{
  Index index        = enzo_block->index();
  const  int level   = index.level();
  // copy face data to FieldFace

  // Pack and send "R" to parent

  int ic3[3];
  index.child(level,&ic3[0],&ic3[1],&ic3[2],min_level_);

  // <COMMON CODE> in restrict_send_() and prolong_send_()

  int if3[3] = {0,0,0};
  int g3[3] = {0,0,0};
  Refresh * refresh = new Refresh;
  refresh->set_prolong(index_prolong_);
  refresh->set_restrict(index_restrict_);
  refresh->add_field(ir_);
  refresh->set_min_face_rank(cello::rank() - 1);

  // copy data from EnzoBlock to array via FieldFace

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, g3, refresh_coarse, refresh);

  refresh->set_restrict(index_restrict_);

  int narray;
  char * array;

  Field field = enzo_block->data()->field();

  field_face->face_to_array(field,&narray,&array);

  delete field_face;

  // Create a FieldMsg for sending data to parent
  // (note: charm messages not deleted on send; are deleted on receive)

  FieldMsg * msg  = new (narray) FieldMsg;

  /// WARNING: double copy

  // Copy FieldFace data to msg

  msg->n = narray;
  memcpy (msg->a, array, narray);
  delete [] array;
  msg->ic3[0] = ic3[0];
  msg->ic3[1] = ic3[1];
  msg->ic3[2] = ic3[2];

  return msg;

}

//----------------------------------------------------------------------

void EnzoSolverMg0::unpack_residual_
(EnzoBlock * enzo_block,FieldMsg * msg) throw()
{
  int if3[3] = {0,0,0};
  int g3[3] = {0,0,0};
  Refresh * refresh = new Refresh;
  refresh->set_prolong(index_prolong_);
  refresh->set_restrict(index_restrict_);
  refresh->add_field(ib_);
  refresh->set_min_face_rank(cello::rank() - 1);

  // copy data from msg to this EnzoBlock

  int * ic3 = msg->ic3;

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, g3, refresh_coarse, refresh);

  refresh->set_restrict(index_restrict_);

  Field field = enzo_block->data()->field();

  char * a = msg->a;
  field_face->array_to_face(a, field);
  delete field_face;

  delete msg;
}

//----------------------------------------------------------------------

FieldMsg * EnzoSolverMg0::pack_correction_
(EnzoBlock * enzo_block, int ic3[3]) throw()
{
  // Pack and send "X" to children

  // <COMMON CODE> in restrict_send_() and prolong_send_()

  int if3[3] = {0,0,0};
  int g3[3];
  cello::field_descr()->ghost_depth(ix_,g3,g3+1,g3+2);
  Refresh * refresh = new Refresh;
  refresh->set_prolong(index_prolong_);
  refresh->set_restrict(index_restrict_);
  refresh->add_field(ix_);
  refresh->set_min_face_rank(cello::rank() - 1);
  // copy data from EnzoBlock to array via FieldFace

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, g3, refresh_fine, refresh);

  Field field = enzo_block->data()->field();
  int narray;
  char * array;
  field_face->face_to_array (field,&narray,&array);

  delete field_face;

  // Create a FieldMsg for sending data to parent
  // (note: charm messages not deleted on send; are deleted on receive)

  FieldMsg * msg  = new (narray) FieldMsg;

  /// WARNING: double copy

  // Copy FieldFace data to msg

  msg->n = narray;
  memcpy (msg->a, array, narray);
  delete [] array;
  msg->ic3[0] = ic3[0];
  msg->ic3[1] = ic3[1];
  msg->ic3[2] = ic3[2];

  //  </COMMON CODE>

  return msg;
}

//----------------------------------------------------------------------

void EnzoSolverMg0::unpack_correction_
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
{
  int if3[3] = {0,0,0};
  int g3[3];
  cello::field_descr()->ghost_depth(ic_,g3,g3+1,g3+2);
  Refresh * refresh = new Refresh;
  refresh->set_prolong(index_prolong_);
  refresh->set_restrict(index_restrict_);
  refresh->add_field(ic_);
  refresh->set_min_face_rank(cello::rank() - 1);

  // copy data from msg to this EnzoBlock

  FieldFace * field_face = enzo_block->create_face
    (if3, msg->ic3, g3, refresh_fine, refresh);

  Field field = enzo_block->data()->field();
  field_face->array_to_face (msg->a, field);

  delete field_face;
  delete msg;

}

//======================================================================

bool EnzoSolverMg0::is_converged_(EnzoBlock * enzo_block) const
{
  return (rr0_ != 0.0 && rr_/rr0_ < res_tol_);
}

//----------------------------------------------------------------------

bool EnzoSolverMg0::is_diverged_(EnzoBlock * enzo_block) const
{
  const int iter = *(((EnzoSolverMg0 *)this)->piter(enzo_block));
  return (iter >= iter_max_);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::end(Block * block)
{
  deallocate_temporary_(block);

  Solver::end_(block);
}
