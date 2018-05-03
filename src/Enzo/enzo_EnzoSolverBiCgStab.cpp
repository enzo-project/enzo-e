// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverBiCgStab.cpp
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-23 16:19:06
/// @brief    Implements the EnzoSolverBiCgStab class

/// The following is Latex code for right-preconditioned BiCgStab,
/// optimized for both memory efficiency and for reduced 'reduction'
/// synchronization (Algorithm 1 in "Analysis and practical use of
/// flexible BICGSTAB", Jie Chen, Lois Curfman McInnes, Hong Zhang,
/// Preprint ANL/MCS-P3039-0912)

///
/// $r_0 = b - A x_0$
/// $\bar{r}_0$ arbitrary
/// $p_0 = r_0$
/// for $j=0,1,\ldots$ until convergence do
///    $\tilde{p}_j = M^{-1}p_j$
///    $\alpha_j = (r_j,\bar{r}_0) / (A\tilde{p}_j,\bar{r}_0)$
///    $s_j = r_j - \alpha_jA\tilde{p}_j$
///    $\tilde{s}_j = M^{-1}s_j$
///    $\omega_j = (A\tilde{s}_j,s_j) / (A\tilde{s}_j,A\tilde{s}_j)$
///    $x_{j+1} = x_j + \alpha_j\tilde{p}_j + \omega_j\tilde{s}_j$
///    $r_{j+1} = s_j - \omega_jA\tilde{s}_j$
///    $\beta_j = (r_{j+1},\bar{r}_0) / (r_j,\bar{r}_0) \cdot \alpha_j/\omega_j$
///    $p_{j+1} = r_{j+1} + \beta_j(p_j - \omega_jA\tilde{p}_j)$
/// end for

/// Below is this algorithm with variable name changes as implemented
/// in EnzoSolverBiCgStab.

/// LINE 01:  R0 = B - A * X_0
/// LINE 02:  P0 = R0
/// LINE 03:  for j=0,1,... until convergence
/// LINE 04:     Y = M \ P
/// LINE 05:     V = A * Y
/// LINE 06:     beta_n = R*R0
/// LINE 07:     alpha = beta_n / (V*R0)
/// LINE 08:     Q = R - alpha * V
/// LINE 09:     X = X + alpha * Y
/// LINE 10:     Y = M \ Q
/// LINE 11:     U = A * Y
/// LINE 12:     omega = (U*Q) / (U*U)
/// LINE 13:     X = X + omega * Y
/// LINE 14:     R = Q - omega * U
/// LINE 15:     beta = (R*R0) / beta_n * (alpha/omega)
/// LINE 16:     P = R + beta * (P - omega * V)
/// LINE 17:  end for

#include "cello.hpp"
#include "charm_simulation.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSolverBiCgStab::EnzoSolverBiCgStab
(FieldDescr* field_descr,
 int monitor_iter, int restart_cycle, int rank,
 int iter_max, double res_tol,
 int min_level, int max_level,
 int index_precon
 ) 
  : Solver(monitor_iter,restart_cycle,min_level,max_level), 
    A_(NULL),
    index_precon_(index_precon),
    rank_(rank),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    rho0_(0), err_(0), err0_(0),err_min_(0), err_max_(0),
    ib_(0), ix_(0), ir_(0), ir0_(0), ip_(0), 
    iy_(0), iv_(0), iq_(0), iu_(0),
    m_(0), mx_(0), my_(0), mz_(0),
    gx_(0), gy_(0), gz_(0),
    iter_(0),
    beta_d_(0), beta_n_(0), beta_(0), 
    omega_d_(0), omega_n_(0), omega_(0), 
    vr0_(0), rr_(0), alpha_(0),
    bs_(0.0),xs_(0.0),c_(0.0),
    ys_(0.0),vs_(0.0),us_(0.0),
    bnorm_(0.0)
{

  ir_ = field_descr->insert_temporary();
  ir0_ = field_descr->insert_temporary();
  ip_ = field_descr->insert_temporary();
  iy_ = field_descr->insert_temporary();
  iv_ = field_descr->insert_temporary();
  iq_ = field_descr->insert_temporary();
  iu_ = field_descr->insert_temporary();
  
  /// Initialize default Refresh (called before entry to compute())
  
  const int ir = add_refresh(4, 0, neighbor_type_(),
			     sync_type_(),
			     enzo_sync_id_solver_bicgstab);
  
  refresh(ir)->add_all_fields();
  
  refresh(ir)->add_field (ir_);
  refresh(ir)->add_field (ir0_);
  refresh(ir)->add_field (ip_);
  refresh(ir)->add_field (iy_);
  refresh(ir)->add_field (iv_);
  refresh(ir)->add_field (iq_);
  refresh(ir)->add_field (iu_);
  
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::apply
( std::shared_ptr<Matrix> A, int ix, int ib, Block * block) throw()
{

  Solver::begin_(block);
  
  A_ = A;
  ix_ = ix;
  ib_ = ib;

  Field field = block->data()->field();

  allocate_temporary_(field,block);

  /// cast input argument to the EnzoBlock associated with this char

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (block);

  /// access the field infromation on this block
  
  field.dimensions (0, &mx_, &my_, &mz_);
  field.ghost_depth(0, &gx_, &gy_, &gz_);

  m_ = mx_*my_*mz_;

  compute_ (enzo_block);
}

//======================================================================

void EnzoSolverBiCgStab::compute_(EnzoBlock* enzo_block) throw() {

  /// initialize BiCgStab iteration counter
  iter_ = 0;

  /// access field container on this block

  Field field = enzo_block->data()->field();

  /// construct RHS B, initialize initial solution X to zero (only on
  /// leaf blocks)

  /// access relevant fields

  enzo_float* X   = (enzo_float*) field.values(ix_);
  enzo_float* R   = (enzo_float*) field.values(ir_);
  enzo_float* R0  = (enzo_float*) field.values(ir0_);
  enzo_float* P   = (enzo_float*) field.values(ip_);
  enzo_float* Y   = (enzo_float*) field.values(iy_);
  enzo_float* V   = (enzo_float*) field.values(iv_);
  enzo_float* Q   = (enzo_float*) field.values(iq_);
  enzo_float* U   = (enzo_float*) field.values(iu_);
  enzo_float* R0_copy  = (enzo_float*) field.values("R0");

  for (int i=0; i<m_; i++) {
    X[i] = 0.0;
    R[i] = 0.0;
    R0[i] = 0.0;
    P[i] = 0.0;
    Y[i] = 0.0;
    V[i] = 0.0;
    Q[i] = 0.0;
    U[i] = 0.0;
    R0_copy[i] = 0.0;
  }
  
  if (is_active_(enzo_block)) {

    if ( reuse_solution_ (enzo_block->cycle()) ) {

      enzo_float* X_copy  = (enzo_float*) field.values("X_copy");

      for (int i=0; i<m_; i++) X[i] = X_copy[i];

      A_->residual (ir_, ib_, ix_, enzo_block);
      
    }
  }

  /// for singular Poisson problems, N(A) is not empty, so project B
  /// into R(A)

  if (A_->is_singular()) {

    /// set reduce[0] = bs_ = SUM(B)
    /// set reduce[1] = xs_ = SUM(X)
    /// set reduce[2] = c_ = COUNT(B)

    long double reduce[3] = {0.0};
    long double count = 0.0;

    if (is_active_(enzo_block)) {

      enzo_float* B = (enzo_float*) field.values(ib_);
      enzo_float* X = (enzo_float*) field.values(ix_);

      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    count++;
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[0] += B[i];
	    reduce[1] += X[i];
	  }
	}
      }
      reduce[2] = count;
    }

    /// initiate callback for r_solver_bicgstab_start_1, reduce sum
    /// and count over all blocks, and continue with
    /// r_solver_bicgstab_start_1()
    
    CkCallback callback(CkIndex_EnzoBlock::r_solver_bicgstab_start_1(NULL), 
			enzo_block->proxy_array());
    
    enzo_block->contribute(3*sizeof(long double), &reduce, 
			   sum_long_double_3_type, callback);

  } else {

    /// for nonsingular systems, call start_2 directly
    this->start_2(enzo_block);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_start_1(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  /// EnzoBlock accumulates global contributions to SUM(B) and COUNT(B)

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());

  long double* data = (long double*) msg->getData();

  solver->set_bs( data[0] );
  solver->set_xs( data[1] );
  solver->set_c ( data[2] );

  delete msg;

  solver->start_2(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::start_2(EnzoBlock* enzo_block) throw() {

  /// access field container on this block

  Field field = enzo_block->data()->field();

  /// update B and initialize temporary vectors (on leaf blocks only)

  long double reduce[2] = {0.0};

  if (is_active_(enzo_block)) {

    /// access relevant fields
    enzo_float* B  = (enzo_float*) field.values(ib_);
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    enzo_float* P  = (enzo_float*) field.values(ip_);
    enzo_float* R  = (enzo_float*) field.values(ir_);
    enzo_float* X  = (enzo_float*) field.values(ix_);

    /// for singular problems, project B into R(A)

    if (A_->is_singular()) {
      enzo_float b_shift = bs_ / c_;
      enzo_float x_shift = xs_ / c_;

      for (int i=0; i<m_; i++) {
	B[i] -= b_shift;
	X[i] -= x_shift;
      }
    }

    // recompute residual given shifted B and X
    A_->residual (ir_, ib_, ix_, enzo_block);
      
    for (int i=0; i<m_; i++) {
      R0[i] = R[i];
      P[i]  = R[i];
    }      

    /// Compute local contributions to beta_n_ = DOT(R, R0)
    /// and B*B for stopping criteria
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[0] += R[i]*R0[i];
	  reduce[1] += B[i]*B[i];
	}
      }
    }
  }
  
  /// initiate callback for r_solver_bicgstab_start_3, compute reductions
  /// over all blocks, and continue with r_solver_bicgstab_start_3()

  CkCallback callback(CkIndex_EnzoBlock::r_solver_bicgstab_start_3(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute(2*sizeof(long double), &reduce, 
			 sum_long_double_2_type, callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_start_3(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  /// EnzoBlock accumulates global contributions to DOT(R, R)

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());

  long double* data = (long double*) msg->getData();

  solver->set_beta_n( data[0] );
  solver->set_bnorm ( data[1] );

  delete msg;

  solver->loop_0(this);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_0(EnzoBlock* enzo_block) throw() {

  /// verify legal floating-point value for preceding reduction result
  
  cello::check(beta_n_,"BCG_beta_n_",__FILE__,__LINE__);

  /// initialize/update current error, store error statistics
  
  const int cycle = enzo_block->cycle();
  bool reuse_x = reuse_solution_(cycle);

  if (iter_ == 0) {
    rho0_ = sqrt(bnorm_); // ||B||
    if (rho0_ == 0.0)  rho0_ = 1.0;
    err_ = sqrt(beta_n_) / rho0_;
    err0_ = err_;
    err_min_ = err_;
    err_max_ = err_;
  } else {
    err_ = sqrt(rr_) / rho0_;
    err_min_ = std::min(err_, err_min_);
    err_max_ = std::max(err_, err_max_);
  }

  const bool is_converged = (err_ < res_tol_);
  const bool is_diverged  = (iter_ >= iter_max_);
  
  /// monitor output solution progress (iteration, residual, etc)

  const bool l_output =
    ( enzo_block->index().is_root() &&
      ( (iter_ == 0) ||
	(is_converged || is_diverged) ||
	(monitor_iter_ && (iter_ % monitor_iter_) == 0 )) );

  if (l_output) {
    if (reuse_x) {
      monitor_output_(enzo_block,iter_,sqrt(bnorm_),
		      err_min_*sqrt(bnorm_),
		      err_*sqrt(bnorm_),
		      err_max_*sqrt(bnorm_));
    } else {
      monitor_output_(enzo_block,iter_,err0_,err_min_,err_,err_max_);
    }
  }

  /// Write final status if done
  if (enzo_block->index().is_root() && (is_converged || is_diverged) ) {
    CkPrintf ("%s DEBUG_SOLVER bicgstab "
	      "final iter = %d rr = %Lg  rho0 = %Lg  rr/rho0 = %Lg\n",
	      enzo_block->name().c_str(),
	      iter_,rr_,rho0_,sqrt(rr_)/rho0_);
  }

  if (is_converged) {

    /// Save copy of X if it will be used as initial guess next cycle

    bool reuse_next_x = reuse_solution_(cycle + 1);

    if (reuse_next_x) {

      Field field = enzo_block->data()->field();

      enzo_float* X  = (enzo_float*) field.values(ix_);
      enzo_float* X_copy  = (enzo_float*) field.values("X_copy");

      for (int i=0; i<m_; i++) X_copy[i] = X[i];
      
    }

    /// Exit the solver with converged status
    this->end(enzo_block, return_converged);

  } else if (is_diverged)  {

    /// Exit the solver with diverged status
    this->end(enzo_block, return_diverged);

    
  } else {

    // else continue with loop_2()

    loop_2(enzo_block);
  
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_2(EnzoBlock* enzo_block) throw() {

  /// access field container on this block

  Field field = enzo_block->data()->field();

  if (index_precon_ >= 0) {

    for (int i=0; i<m_; i++) Y[i] = 0.0;

    /// Access the preconditioner for this solver, if any

    Simulation * simulation = proxy_simulation.ckLocalBranch();
    Solver * precon = simulation->problem()->solver(index_precon_);

    /// Apply the preconditioner, then return to
    /// p_solver_bicgstab_loop_2()
    
    precon->set_sync_id (8);
    precon->set_callback(CkIndex_EnzoBlock::p_solver_bicgstab_loop_2());

    /// LINE 04: Y = M \ P
    
    precon->apply(A_,iy_,ip_,enzo_block);
    
  } else { // no preconditioner

    enzo_float * Y = (enzo_float*) field.values(iy_);
    enzo_float * P = (enzo_float*) field.values(ip_);
    
    /// LINE 04: Y = M \ P  [ M = I ]
    
    for (int i=0; i<m_; i++) Y[i] = P[i];

    loop_25(enzo_block);

  }

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_2() {

  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  solver->loop_25(enzo_block);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_25 (EnzoBlock * enzo_block) throw() {
  
  // Refresh field faces then call p_solver_bicgstab_loop_25()
  
  Refresh refresh (4,0,neighbor_type_(), sync_type_(),
		   enzo_sync_id_solver_bicgstab_loop_25);
  
  refresh.set_active(is_active_(enzo_block));
  
  refresh.add_all_fields();

  refresh.add_field (ir_);
  refresh.add_field (ir0_);
  refresh.add_field (ip_);
  refresh.add_field (iy_);
  refresh.add_field (iv_);
  refresh.add_field (iq_);
  refresh.add_field (iu_);

  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_bicgstab_loop_3(),&refresh);

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_3() {

  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  solver->loop_4(enzo_block);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
  
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_4(EnzoBlock* enzo_block) throw() {

  /// access field container on this block

  Field field = enzo_block->data()->field();

  /// V = MATVEC(A,Y)
  
  if (is_active_(enzo_block)) {

    /// LINE 05: V = A * Y
    
    A_->matvec(iv_, iy_, enzo_block);     

  }

  /// compute local contributions to vr0_ = DOT(V, R0)
  
  long double reduce[3] = {0.0};
  
  if (is_active_(enzo_block)) {
    
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    enzo_float* V  = (enzo_float*) field.values(iv_);

    /// LINE 07 [part]  vr0_ = V*R0
    
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[0] += V[i]*R0[i];
	}
      }
    }
    
  }

  /// for singular Poisson problems need all vectors in R(A), so
  /// project both Y and V into R(A)

  if (A_->is_singular()) {

    if (is_active_(enzo_block)) {

      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* V = (enzo_float*) field.values(iv_);

      /// ys_ = sum (Y[i])
      /// vs_ = sum (V[i])
      
      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[1] += Y[i];
	    reduce[2] += V[i];
	  }
	}
      }
    }
  }

  /// contribute to global sums over blocks, and return
  /// r_solver_bicgstab_loop_5()

  CkCallback callback(CkIndex_EnzoBlock::r_solver_bicgstab_loop_5(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute(3*sizeof(long double), &reduce, 
			 sum_long_double_3_type, callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_5(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());

  long double* data = (long double*) msg->getData();
  
  solver->set_vr0( data[0] );
  solver->set_ys(  data[1] );
  solver->set_vs(  data[2] );
  
  delete msg;

  solver->loop_6(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_6(EnzoBlock* enzo_block) throw() {

  /// access field container on this block

  Field field = enzo_block->data()->field();

  /// for singular problems, project Y and V into R(A)
  
  if (is_active_(enzo_block) && A_->is_singular()) {
    
    enzo_float* Y = (enzo_float*) field.values(iy_);
    enzo_float* V = (enzo_float*) field.values(iv_);

    enzo_float y_shift = ys_ / c_;
    enzo_float v_shift = vs_ / c_;

    for (int i=0; i<m_; i++) {
      Y[i] -= y_shift;
      V[i] -= v_shift;
    }
  }

  /// compute alpha factor in BiCgStab algorithm (all blocks)

  alpha_ = beta_n_ / vr0_;

  /// update vectors (on leaf blocks only)

  if (is_active_(enzo_block)) {

    /// access relevant fields
    enzo_float* Q = (enzo_float*) field.values(iq_);
    enzo_float* R = (enzo_float*) field.values(ir_);
    enzo_float* V = (enzo_float*) field.values(iv_);
    enzo_float* X = (enzo_float*) field.values(ix_);
    enzo_float* Y = (enzo_float*) field.values(iy_);

    /// LINE 08: Q = R - alpha * V
    /// LINE 09: X = X + alpha * Y
    
    for (int i=0; i<m_; i++) {
      Q[i] = R[i] - alpha_*V[i];
      X[i] = X[i] + alpha_*Y[i];
    }
  }

  /// refresh Q with callback to p_solver_bicgstab_loop_7 to handle re-entry

  loop_8(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_8(EnzoBlock* enzo_block) throw() {

  /// access field container on this block

  Field field = enzo_block->data()->field();

  if (index_precon_ >= 0) {

    enzo_float* Y = (enzo_float*) field.values(iy_);

    for (int i=0; i<m_; i++) Y[i] = 0.0;

    /// Access the preconditioner for this solver, if any

    Simulation * simulation = proxy_simulation.ckLocalBranch();
    Solver * precon = simulation->problem()->solver(index_precon_);

    /// Apply the preconditioner, then return to
    /// p_solver_bicgstab_loop_8()
    
    precon->set_sync_id (10);
    precon->set_callback(CkIndex_EnzoBlock::p_solver_bicgstab_loop_8());

    /// LINE 10: Y = M \ Q
    
    precon->apply(A_,iy_,iq_,enzo_block);
    
  } else {

    enzo_float * Y = (enzo_float*) field.values(iy_);
    enzo_float * Q = (enzo_float*) field.values(iq_);

    /// LINE 10: Y = M \ Q  [ M = I ]

    for (int i=0; i<m_; i++)  Y[i] = Q[i];

    loop_85(enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_8() {

  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  
  solver->loop_85(enzo_block);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_85 (EnzoBlock * enzo_block) throw() {
  
  // Refresh field faces then call p_solver_bicgstab_loop_85()

  Refresh refresh (4,0,neighbor_type_(), sync_type_(),
		   enzo_sync_id_solver_bicgstab_loop_85);
  
  refresh.set_active(is_active_(enzo_block));
  refresh.add_all_fields();

  refresh.add_field (ir_);
  refresh.add_field (ir0_);
  refresh.add_field (ip_);
  refresh.add_field (iy_);
  refresh.add_field (iv_);
  refresh.add_field (iq_);
  refresh.add_field (iu_);
  
  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_bicgstab_loop_9(),&refresh);

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_9() {

  performance_start_(perf_compute,__FILE__,__LINE__);
  
  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  solver->loop_10(enzo_block);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
  
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_10(EnzoBlock* enzo_block) throw() {
  
  /// access field container on this block

  Field field = enzo_block->data()->field();

  if (is_active_(enzo_block)) {

    /// LINE 11:     U = A * Y
    
    A_->matvec(iu_, iy_, enzo_block);     /// apply matrix to local block

  }

  long double reduce[4] = {0.0};
  
  if (is_active_(enzo_block)) {
    
    enzo_float* U  = (enzo_float*) field.values(iu_);
    enzo_float* Q  = (enzo_float*) field.values(iq_);
    
    /// omega_d_ = DOT(U, U)
    /// omega_n_ = DOT(U, Q)
  
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[0] += U[i]*U[i];
	  reduce[1] += U[i]*Q[i];
	}
      }
    }
  
    /// for singular Poisson problems, project both Y and U into R(A)

    if (A_->is_singular()) {

      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* U = (enzo_float*) field.values(iu_);

      /// ys_ = SUM(Y)
      /// us_ = SUM(U)

      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[2] += Y[i];
	    reduce[3] += U[i];
	  }
	}
      }
    }
  }
  
  /// compute sums over Blocks and continue with r_solver_bicgstab_loop_11()

  CkCallback callback(CkIndex_EnzoBlock::r_solver_bicgstab_loop_11(NULL), 
		      enzo_block->proxy_array());
  
  enzo_block->contribute(4*sizeof(long double), &reduce, 
			 sum_long_double_4_type, callback);
    
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_11(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());
  
  long double* data = (long double*) msg->getData();
  
  solver->set_omega_d( data[0] );
  solver->set_omega_n( data[1] );
  solver->set_ys( data[2] );
  solver->set_us( data[3] );
  
  delete msg;

  solver->loop_12(this);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_12(EnzoBlock* enzo_block) throw() {

  /// verify legal floating-point values for preceding reduction results

  cello::check(omega_d_,"BCG_omega_d_",__FILE__,__LINE__);
  cello::check(omega_n_,"BCG_omega_n_",__FILE__,__LINE__);

  /// access field container on this block

  Field field = enzo_block->data()->field();

  /// for singular problems, update omega_d_ and project Y and U into R(A)

  if (A_->is_singular()) {

    omega_d_ = omega_d_ - us_*us_/c_;

    if (is_active_(enzo_block)) {
      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* U = (enzo_float*) field.values(iu_);
      enzo_float y_shift = ys_ / c_;
      enzo_float u_shift = us_ / c_;
      for (int i=0; i<m_; i++) {
	Y[i] -= y_shift;
	U[i] -= u_shift;
      }
    }
  }

  /// avoid division by 0.0
  
  if (omega_d_ == 0.0)  omega_d_ = 1.0;
  
  /// LINE 12:     omega = (U*Q) / (U*U)
  
  omega_ = omega_n_ / omega_d_;

  /// check for breakdown in BiCgStab
  
  if ( omega_ == 0.0 ) {
    WARNING ("EnzoSolverBiCgStab::loop12()",
	     "Solver error: omega_ == 0");
    this->end(enzo_block, return_error);
  }

  /// update vectors on leaf blocks
  
  if (is_active_(enzo_block)) {

    enzo_float* X = (enzo_float*) field.values(ix_);
    enzo_float* Y = (enzo_float*) field.values(iy_);
    enzo_float* R = (enzo_float*) field.values(ir_);
    enzo_float* Q = (enzo_float*) field.values(iq_);
    enzo_float* U = (enzo_float*) field.values(iu_);
    
    /// LINE 13:     X = X + omega * Y
    /// LINE 14:     R = Q - omega * U

    for (int i=0; i<m_; i++) {
      X[i] = X[i] + omega_*Y[i];
      R[i] = Q[i] - omega_*U[i];
    }
  }

  
  /// Update previous beta value (beta_d_) to current value (beta_n_)
  
  beta_d_ = beta_n_;
  
  /// rr_     = DOT(R, R)
  /// beta_n_ = DOT(R, R0)
  
  long double reduce[2] = {0.0};
  
  if (is_active_(enzo_block)) {
    
    enzo_float* R  = (enzo_float*) field.values(ir_);
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[0] += R[i]*R[i];
	  reduce[1] += R[i]*R0[i];
	}
      }
    }
  }

  /// sum over blocks and continue with r_solver_bicgstab_loop_13()
  
  CkCallback callback(CkIndex_EnzoBlock::r_solver_bicgstab_loop_13(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute(2*sizeof(long double), &reduce, 
			 sum_long_double_2_type, callback);
    
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_13(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());
  
  long double* data = (long double*) msg->getData();
  
  solver->set_rr(     data[0] );
  solver->set_beta_n( data[1] );
  
  delete msg;

  solver->loop_14(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_14(EnzoBlock* enzo_block) throw() {

  /// verify legal floating-point values for preceding reduction results
  
  cello::check(rr_,    "BCG_rr_",    __FILE__,__LINE__);
  cello::check(beta_n_,"BCG_beta_n_",__FILE__,__LINE__);

  /// access field container on this block

  Field field = enzo_block->data()->field();

  /// check for breakdown in BiCgStab
  
  if (beta_n_ == 0.0) {
    WARNING ("EnzoSolverBiCgStab::loop14()",
	     "Solver error: beta_n_ == 0");
    this->end(enzo_block, return_error);
  }
  

  /// LINE 15: beta = (R*R0) / beta_n * (alpha/omega)
  
  beta_ = (beta_n_/beta_d_)*(alpha_/omega_);

  /// update direction vector
  
  if (is_active_(enzo_block)) {

    enzo_float* P = (enzo_float*) field.values(ip_);
    enzo_float* R = (enzo_float*) field.values(ir_);
    enzo_float* V = (enzo_float*) field.values(iv_);

    /// LINE 16:     P = R + beta * (P - omega * V)

    for (int i=0; i<m_; i++) {
      P[i] = R[i] + beta_*(P[i] - omega_*V[i]);
    }
  }

  /// contribute to global iteration counter and continue with
  /// r_solver_bicgstab_loop_15()
  
  int iter = iter_ + 1;

  CkCallback callback(CkIndex_EnzoBlock::r_solver_bicgstab_loop_15(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute(sizeof(int), &iter, 
			 CkReduction::max_int, callback);

}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_15(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverBiCgStab* solver = 
    static_cast<EnzoSolverBiCgStab*> (this->solver());
  
  solver->set_iter ( ((int*)msg->getData())[0] );
  
  delete msg;

  solver->loop_0(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::end (EnzoBlock* enzo_block, int retval) throw () {

  Field field = enzo_block->data()->field();

  deallocate_temporary_(field,enzo_block);
  
  Solver::end_(enzo_block);
  
  CkCallback(callback_,
	     CkArrayIndexIndex(enzo_block->index()),
	     enzo_block->proxy_array()).send();

}

//======================================================================
