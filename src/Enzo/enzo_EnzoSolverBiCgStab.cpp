// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverBiCgStab.cpp
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-23 16:19:06
/// @brief    Implements the EnzoSolverBiCgStab class

/// Below is this algorithm with variable name changes as implemented
/// in EnzoSolverBiCgStab.  This is based on Algorithm 1 in "Analysis
/// and practical use of flexible BICGSTAB", Jie Chen, Lois Curfman
/// McInnes, Hong Zhang, Preprint ANL/MCS-P3039-0912)

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

// #define TRACE_SOLVER_BCG

// #define DEBUG_FIELD
// #define DEBUG_GHOST
// #define DEBUG_SCALAR

// #define DEBUG_WRITE
// #define DEBUG_READ

#ifdef NEW_CONTRIBUTE
static const int TEMP_INC = 1;
#else
static const int TEMP_INC = 0;
#endif

#ifdef DEBUG_SCALAR
#   define TRACE_SCALAR(BLOCK,NAME,SCALAR)				\
  CkPrintf ("%s:%d %s TRACE_SCALAR %s = %20.15lg\n",			\
	    __FILE__,__LINE__,BLOCK->name().c_str(),NAME,SCALAR);
  #else
#   define TRACE_SCALAR(BLOCK,NAME,SCALAR) /* ... */
#endif

#ifdef DEBUG_GHOST
#   define TRACE_GHOST(BLOCK,ID,NAME)					\
  {									\
    Field field = block->data()->field();				\
    enzo_float* X      = (enzo_float*) field.values(ID);		\
    double sum_a[4],sum_aa[4],sum_abs[4];				\
    for (int i=0; i<4; i++) sum_a[i] = 0.0;				\
    for (int i=0; i<4; i++) sum_aa[i] = 0.0;				\
    for (int i=0; i<4; i++) sum_abs[i] = 0.0;				\
    int ig = 1;								\
    int ix,iy,iz;							\
    int iz0 = gz_;							\
    int iz1 = gz_-ig;							\
    for (int iy=gy_; iy<my_-gy_; iy++) {				\
      for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	int i0=ix+mx_*(iy+my_*iz0);					\
	int i1=ix+mx_*(iy+my_*iz1);					\
	enzo_float v0 = X[i0];						\
	enzo_float v1 = X[i1];						\
	sum_a[ig-1]   += (v0-v1);					\
	sum_aa[ig-1]  += (v0-v1)*(v0-v1);				\
	sum_abs[ig-1] += std::abs(v0-v1);				\
      }									\
    }									\
    iz0 = mz_-gz_-1;							\
    iz1 = mz_-gz_-1+ig;							\
    for (int iy=gy_; iy<my_-gy_; iy++) {				\
      for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	int i0=ix+mx_*(iy+my_*iz0);					\
	int i1=ix+mx_*(iy+my_*iz1);					\
	enzo_float v0 = X[i0];						\
	enzo_float v1 = X[i1];						\
	sum_a[ig-1]   += (v0-v1);					\
	sum_aa[ig-1]  += (v0-v1)*(v0-v1);				\
	sum_abs[ig-1] += std::abs(v0-v1);				\
      }									\
    }									\
    CkPrintf ("%s:%d %s TRACE_GHOST %s layer %d sum(A) sum(A*A) sum(abs(A)) " \
	      " %20.15lg %20.15lg %20.15lg\n"				\
	      ,__FILE__,__LINE__,BLOCK->name().c_str(),NAME,ig,sum_a[0],sum_aa[0],sum_abs[0]); \
  }
#else
#   define TRACE_GHOST(BLOCK,ID,NAME) /* ... */
#endif


#ifdef DEBUG_FIELD
#   define COPY_FIELD(BLOCK,ID,COPY)					\
  {									\
    Field field = block->data()->field();				\
    enzo_float* X      = (enzo_float*) field.values(ID);		\
    enzo_float* X_bcg  = (enzo_float*) field.values(COPY);		\
    if (X_bcg) for (int i=0; i<m_; i++)  X_bcg[i] = X[i];		\
    double sum_a=0.0,sum_abs=0.0;					\
    for (int iz=gz_; iz<mz_-gz_; iz++) {				\
      for (int iy=gy_; iy<my_-gy_; iy++) {				\
	for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	  int i=ix+mx_*(iy+my_*iz);					\
	  sum_a+=X[i];							\
	  sum_abs+=std::abs(X[i]);					\
	}								\
      }									\
    }									\
    CkPrintf ("%s:%d %s %s COPY_FIELD %d %s shift %20.15lg %20.15lg\n"	\
	      ,__FILE__,__LINE__,BLOCK->name().c_str(),name().c_str(),ID,COPY,sum_a, sum_abs); \
    TRACE_GHOST(BLOCK,ID,COPY);						\
  }
#else
#   define COPY_FIELD(BLOCK,ID,COPY) /* ... */
#endif

  
//----------------------------------------------------------------------

EnzoSolverBiCgStab::EnzoSolverBiCgStab
(std::string name,
 std::string field_x, std::string field_b,
 int monitor_iter, int restart_cycle,
 int solve_type,
 int min_level, int max_level,
 int iter_max, double res_tol,
 int index_precon
 ) 
  : Solver(name,
	   field_x,
	   field_b,
	   monitor_iter,
	   restart_cycle,
	   solve_type,
	   min_level,
	   max_level),
    alpha_(0), beta_n_(0), beta_d_(0),   omega_(0),
    rr_(0), r0s_(0.0), c_(0.0), bnorm_(0.0),
    rho0_(0), err_(0), err0_(0), err_min_(0), err_max_(0),
    res_tol_(res_tol),
    A_(NULL),
    index_precon_(index_precon),
    iter_max_(iter_max), 
    ir_(0), ir0_(0), ip_(0), 
    iy_(0), iv_(0), iq_(0), iu_(0),
    nx_(0), ny_(0), nz_(0),
    m_(0), mx_(0), my_(0), mz_(0),
    gx_(0), gy_(0), gz_(0),
    iter_(0)
{
  FieldDescr * field_descr = cello::field_descr();

  CkPrintf ("%s:%d DEBUG_NEW_CONTRIBUTE solver %s solve_type %d\n",
	    __FILE__,__LINE__,this->name().c_str(),solve_type);
  if (solve_type == solve_tree) {
    ScalarDescr * scalar_descr = cello::scalar_descr_long_double();
    is_alpha_ =  scalar_descr->new_value("solver_bicgstab_alpha");
    is_beta_n_ = scalar_descr->new_value("solver_bicgstab_beta_n");
    is_beta_d_ = scalar_descr->new_value("solver_bicgstab_beta_d");
    is_omega_ =  scalar_descr->new_value("solver_bicgstab_omega");
    is_rr_ =     scalar_descr->new_value("solver_bicgstab_rr");
    is_r0s_ =    scalar_descr->new_value("solver_bicgstab_r0s");
    is_c_ =      scalar_descr->new_value("solver_bicgstab_c");
    is_bnorm_ =  scalar_descr->new_value("solver_bicgstab_bnorm");
    is_vr0_ =    scalar_descr->new_value("solver_bicgstab_vr0");
    is_ys_ =     scalar_descr->new_value("solver_bicgstab_ys");
    is_vs_ =     scalar_descr->new_value("solver_bicgstab_vs");
    is_omega_d_ =scalar_descr->new_value("solver_bicgstab_omega_d");
    is_omega_n_ =scalar_descr->new_value("solver_bicgstab_omega_n");
    is_us_ =     scalar_descr->new_value("solver_bicgstab_us");
    is_qs_ =     scalar_descr->new_value("solver_bicgstab_qs");
  }
  
  ir_ = field_descr->insert_temporary();
  ir0_ = field_descr->insert_temporary();
  ip_ = field_descr->insert_temporary();
  iy_ = field_descr->insert_temporary();
  iv_ = field_descr->insert_temporary();
  iq_ = field_descr->insert_temporary();
  iu_ = field_descr->insert_temporary();
  
  /// Initialize default Refresh (called before entry to compute())

  const int min_face_rank = cello::rank() - 1;
  
  const int ir = add_refresh(4, min_face_rank, neighbor_type_(),
			     sync_type_(),
			     enzo_sync_id_solver_bicgstab);

  refresh(ir)->add_field (field_x);
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
( std::shared_ptr<Matrix> A, Block * block) throw()
{

  Solver::begin_(block);
  
  A_ = A;

  Field field = block->data()->field();

  allocate_temporary_(block);

  /// cast input argument to the EnzoBlock associated with this char

  EnzoBlock* enzo_block = enzo::block(block);

  /// access the field infromation on this block
  
  field.dimensions (0, &mx_, &my_, &mz_);
  field.ghost_depth(0, &gx_, &gy_, &gz_);

  m_ = mx_*my_*mz_;

  compute_ (enzo_block);
}

//======================================================================

void EnzoSolverBiCgStab::compute_(EnzoBlock* block) throw() {

  /// initialize BiCgStab iteration counter
  iter_ = 0;

  /// access field container on this block

  Field field = block->data()->field();

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

  for (int i=0; i<m_; i++) {
    X[i] = R[i] = R0[i] = P[i] = 0.0;
    Y[i] = V[i] = Q[i] =  U[i] = 0.0;
  }
  
  if (is_finest_(block)) {

    if ( reuse_solution_ (block->cycle()) ) {

      enzo_float* X_copy  = (enzo_float*) field.values("X_copy");
      for (int i=0; i<m_; i++) X[i] = X_copy[i];

      A_->residual (ir_, ib_, ix_, block);
      
    }
  }

  /// for singular Poisson problems, N(A) is not empty, so project B
  /// into R(A)

  if (A_->is_singular()) {

    /// set reduce[0] = bs_ = SUM(B)
    /// set reduce[1] = xs_ = SUM(X)
    /// set reduce[2] = c_ = COUNT(B)

    std::vector<long double> reduce;
    reduce.resize(3+TEMP_INC);
    reduce.clear();
    if (TEMP_INC) reduce[0] = 3;
    
    long double count = 0.0;

    if (is_finest_(block)) {

      enzo_float* B = (enzo_float*) field.values(ib_);
      enzo_float* X = (enzo_float*) field.values(ix_);
 
      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    count++;
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[1+TEMP_INC] += B[i];
	    reduce[2+TEMP_INC] += X[i];
	  }
	}
      }
      reduce[0+TEMP_INC] = count;
    }

    /// initiate callback for r_solver_bicgstab_start_1, reduce sum
    /// and count over all blocks, and continue with
    /// r_solver_bicgstab_start_1()
    
    std::vector<int> is_array;

    if (solve_type_ == solve_tree) {
      is_array.resize(3);
      is_array[0] = is_c_;
      is_array[1] = is_bs_;
      is_array[2] = is_xs_;
    }

    CkCallback callback = CkCallback
      (CkIndex_EnzoBlock::r_solver_bicgstab_start_1(NULL),
       block->proxy_array());

    inner_product_(block,3,&reduce[0],is_array,callback);

  } else {

    this->start_2(block,NULL);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_start_1(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->start_2(this,msg);

  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::start_2(EnzoBlock* block,
				 CkReductionMsg *msg) throw() {

  double bs = 0.0;
  double xs = 0.0;
  c_ = 1.0;

  if (msg != NULL) {

    long double* data = (long double*) msg->getData();

    if (solve_type_ != solve_tree) {
      c_ = data[0+TEMP_INC];
      bs = data[1+TEMP_INC];
      xs = data[2+TEMP_INC];
    } else {
      c_ = *s_c_(block);
      bs = *s_bs_(block);
      xs = *s_xs_(block);
    }

  }

  delete msg;
  
  /// access field container on this block

  Field field = block->data()->field();

  /// update B and initialize temporary vectors (on leaf blocks only)

  std::vector<long double> reduce;
  reduce.resize(3+TEMP_INC);
  reduce.clear();
  if (TEMP_INC) reduce[0] = 3;

  if (is_finest_(block)) {

    /// access relevant fields
    enzo_float* B  = (enzo_float*) field.values(ib_);
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    enzo_float* P  = (enzo_float*) field.values(ip_);
    enzo_float* R  = (enzo_float*) field.values(ir_);
    enzo_float* X  = (enzo_float*) field.values(ix_);

    /// for singular problems, project B into R(A)

    if (A_->is_singular()) {
      enzo_float b_shift = bs / c_;
      enzo_float x_shift = xs / c_;

      for (int i=0; i<m_; i++) {
	B[i] -= b_shift;
	X[i] -= x_shift;
      }
    }

    // recompute residual given shifted B and X
    A_->residual (ir_, ib_, ix_, block);

    /// LINE 01:  R0 = B - A * X_0
    /// LINE 02:  P0 = R0
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
	  reduce[0+TEMP_INC] += R[i]*R0[i];
	  reduce[1+TEMP_INC] += B[i]*B[i];
	  reduce[2+TEMP_INC] += R[i];
	}
      }
    }
  }
  
  /// initiate callback for r_solver_bicgstab_start_3, compute reductions
  /// over all blocks, and continue with r_solver_bicgstab_start_3()

  std::vector<int> is_array;
  if (solve_type_ == solve_tree) {
    is_array.resize(3);
    is_array[0] = is_beta_n_;
    is_array[1] = is_bnorm_;
    is_array[2] = is_r0s_;
  }

  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_solver_bicgstab_start_3(NULL), 
     block->proxy_array());

  inner_product_(block,3,&reduce[0],is_array,callback);

}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_start_3(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_0a(this,msg);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_0a(EnzoBlock* block,
				 CkReductionMsg * msg) throw() {

  long double* data = (long double*) msg->getData();

  beta_n_ = data[0+TEMP_INC];
  bnorm_  = data[1+TEMP_INC];
  r0s_    = data[2+TEMP_INC];
  TRACE_SCALAR(block,"beta_n_",beta_n_);
  TRACE_SCALAR(block,"bnorm_",bnorm_);
  TRACE_SCALAR(block,"r0s_",r0s_);
  delete msg;
  loop_0(block);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_0b(EnzoBlock* block,
				 CkReductionMsg * msg) throw() {

  iter_ = ((int*)msg->getData())[0];
  delete msg;
  loop_0(block);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_0(EnzoBlock* block) throw() {

/// verify legal floating-point value for preceding reduction result
  
  cello::check(beta_n_,"BCG_beta_n_",__FILE__,__LINE__);

  /// initialize/update current error, store error statistics
  
  const int cycle = block->cycle();
  bool reuse_x = reuse_solution_(cycle);

  if (iter_ == 0) {
    if (is_finest_(block)) {
      if (A_->is_singular()) {
	enzo_float shift = r0s_/c_;
	
	Field field = block->data()->field();
	enzo_float* R0 = (enzo_float*) field.values(ir0_);
	enzo_float* P  = (enzo_float*) field.values(ip_);
	enzo_float* R  = (enzo_float*) field.values(ir_);
	for (int i=0; i<m_; i++) {
	  R0[i] -= shift;
	  R[i]  -= shift;
	  P[i]  -= shift;
	}
	beta_n_ -= r0s_*r0s_/c_;
	TRACE_SCALAR(block,"beta_n_",beta_n_);
      }
    }
    rho0_ = sqrt(bnorm_); // ||B||
    TRACE_SCALAR(block,"rho0_",rho0_);
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
  TRACE_SCALAR(block,"err_",err_);


  const bool is_converged = (err_ < res_tol_);
  const bool is_diverged  = (iter_ >= iter_max_);

  /// monitor output solution progress (iteration, residual, etc)

  const bool l_output =
    ( block->index().is_root() &&
      ( (iter_ == 0) ||
	(is_converged || is_diverged) ||
	(monitor_iter_ && (iter_ % monitor_iter_) == 0 )) );

  if (l_output) {
    if (reuse_x) {
      monitor_output_(block,iter_,sqrt(bnorm_),
		      err_min_*sqrt(bnorm_),
		      err_*sqrt(bnorm_),
		      err_max_*sqrt(bnorm_));
    } else {
      monitor_output_(block,iter_,err0_,err_min_,err_,err_max_);
    }
  }

  /// Write final status if done
  if (block->index().is_root() && (is_converged || is_diverged) ) {
    CkPrintf ("%s DEBUG_SOLVER bicgstab "
	      "final iter = %d rr = %lg  rho0 = %lg  rr/rho0 = %lg\n",
	      block->name().c_str(),
	      iter_,rr_,rho0_,sqrt(rr_)/rho0_);
  }

  if (is_converged) {

    /// Save copy of X if it will be used as initial guess next cycle

    bool reuse_next_x = reuse_solution_(cycle + 1);

    if (reuse_next_x) {

      Field field = block->data()->field();

      enzo_float* X       = (enzo_float*) field.values(ix_);
      enzo_float* X_copy  = (enzo_float*) field.values("X_copy");

      for (int i=0; i<m_; i++) X_copy[i] = X[i];
      
    }

    this->end(block, return_converged);

  } else if (is_diverged)  {

    this->end(block, return_diverged);

    
  } else {

    loop_2(block);
  
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_2(EnzoBlock* block) throw() {

  /// access field container on this block

  Field field = block->data()->field();

  if (index_precon_ >= 0) {

    enzo_float * Y = (enzo_float*) field.values(iy_);

    for (int i=0; i<m_; i++) Y[i] = 0.0;

    /// Access the preconditioner for this solver, if any

    Solver * precon = cello::solver(index_precon_);

    /// Apply the preconditioner, then return to
    /// p_solver_bicgstab_loop_2()

    precon->set_sync_id (8);
    precon->set_callback(CkIndex_EnzoBlock::p_solver_bicgstab_loop_2());

    /// LINE 04: Y = M \ P
#ifdef TRACE_SOLVER_BCG    
    CkPrintf ("%s %s:%d TRACE_SOLVER_BCG calling preconditioner\n",
	      block->name().c_str(),__FILE__,__LINE__);
#endif

#ifdef DEBUG_READ
    char buffer[40];
    sprintf (buffer,"Q.%s",block->name().c_str());
    enzo_float * P = (enzo_float*) field.values(ip_);
    FILE * fp = fopen(buffer,"r");
    for (int i=0; i<m_; i++) fscanf (fp,"%g",P+i);
    fclose(fp);
#endif    
    
    
    COPY_FIELD(block,ip_,"P0_bcg");

    precon->set_field_x(iy_);
    precon->set_field_b(ip_);
    precon->apply(A_,block);
    COPY_FIELD(block,iy_,"Y0_bcg");
    
  } else { // no preconditioner

    enzo_float * Y = (enzo_float*) field.values(iy_);
    enzo_float * P = (enzo_float*) field.values(ip_);
    
    /// LINE 04: Y = M \ P  [ M = I ]
    for (int i=0; i<m_; i++) Y[i] = P[i];

    loop_25(block);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_2() {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_25(this);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_25 (EnzoBlock * block) throw() {
  
  // Refresh field faces then call p_solver_bicgstab_loop_25()
  
  const int min_face_rank = cello::rank() - 1;
  Refresh refresh (4,min_face_rank,neighbor_type_(), sync_type_(),
		   enzo_sync_id_solver_bicgstab_loop_25);

  refresh.set_active(is_finest_(block));
  //  refresh.add_all_fields();

  refresh.add_field (ix_);

  refresh.add_field (ir_);
  refresh.add_field (ir0_);
  refresh.add_field (ip_);
  refresh.add_field (iy_);
  refresh.add_field (iv_);
  refresh.add_field (iq_);
  refresh.add_field (iu_);

  block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_bicgstab_loop_3(),&refresh);

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_3() {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_4(this);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
  
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_4(EnzoBlock* block) throw() {

  /// access field container on this block

  Field field = block->data()->field();

  /// V = MATVEC(A,Y)

  COPY_FIELD(block,ip_,"P1_bcg");
  COPY_FIELD(block,iy_,"Y1_bcg");
  
  if (is_finest_(block)) {

    /// LINE 05: V = A * Y
    
    A_->matvec(iv_, iy_, block);

  }

  COPY_FIELD(block,iv_,"V1_bcg");
  /// compute local contributions to vr0_ = DOT(V, R0)
  
  std::vector<long double> reduce;
  reduce.resize(3+TEMP_INC);
  reduce.clear();
  if (TEMP_INC) reduce[0] = 3;
  
  if (is_finest_(block)) {
    
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    enzo_float* V  = (enzo_float*) field.values(iv_);

    /// LINE 07 [part]  vr0_ = V*R0
    
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[0+TEMP_INC] += V[i]*R0[i];
	}
      }
    }
    
    /// for singular Poisson problems need all vectors in R(A), so
    /// project both Y and V into R(A)

    if (A_->is_singular()) {

      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* V = (enzo_float*) field.values(iv_);

      /// ys_ = sum (Y[i])
      /// vs_ = sum (V[i])
      
      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[1+TEMP_INC] += Y[i];
	    reduce[2+TEMP_INC] += V[i];
	  }
	}
      }
    }
  }

  /// contribute to global sums over blocks, and return
  /// r_solver_bicgstab_loop_5()

  std::vector<int> is_array;
  if (solve_type_ == solve_tree) {
    is_array.resize(3);
    is_array[0] = is_vr0_;
    is_array[1] = is_ys_;
    is_array[2] = is_vs_;
  }

  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_solver_bicgstab_loop_5(NULL), 
     block->proxy_array());

  inner_product_(block,3,&reduce[0],is_array,callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_5(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_6(this,msg);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_6(EnzoBlock* block,
				CkReductionMsg * msg) throw() {

  long double* data = (long double*) msg->getData();
  
  double vr0 = data[0+TEMP_INC];
  double ys  = data[1+TEMP_INC];
  double vs  = data[2+TEMP_INC];
  
  TRACE_SCALAR(block,"vr0",vr0);
  TRACE_SCALAR(block,"ys",ys);
  TRACE_SCALAR(block,"vs",vs);

  COPY_FIELD(block,iv_,"V1_bcg");
  delete msg;

  /// access field container on this block

  Field field = block->data()->field();

  /// for singular problems, project Y and V into R(A)
  
  if (is_finest_(block) && A_->is_singular()) {
    
    enzo_float* Y = (enzo_float*) field.values(iy_);
    enzo_float* V = (enzo_float*) field.values(iv_);

    enzo_float y_shift = ys / c_;
    enzo_float v_shift = vs / c_;

    for (int i=0; i<m_; i++) {
      Y[i] -= y_shift;
      V[i] -= v_shift;
    }
    COPY_FIELD(block,iy_,"Y_shift");
    COPY_FIELD(block,iv_,"V_shift");

  }

  /// compute alpha factor in BiCgStab algorithm (all blocks)

  /// LINE 07:     alpha = beta_n / (V*R0)
  alpha_ = beta_n_ / vr0;

  TRACE_SCALAR(block,"alpha",alpha_);

  /// update vectors (on leaf blocks only)

  if (is_finest_(block)) {

    /// access relevant fields
    enzo_float* Q = (enzo_float*) field.values(iq_);
    enzo_float* R = (enzo_float*) field.values(ir_);
    enzo_float* V = (enzo_float*) field.values(iv_);
    enzo_float* X = (enzo_float*) field.values(ix_);
    enzo_float* Y = (enzo_float*) field.values(iy_);

    /// LINE 08: Q = R - alpha * V
    /// LINE 09: X = X + alpha * Y
    
    COPY_FIELD(block,ir_,"R1_bcg");
    for (int i=0; i<m_; i++) {
      Q[i] = R[i] - alpha_*V[i];
      X[i] = X[i] + alpha_*Y[i];
    }
  }
  COPY_FIELD(block,iq_,"Q");
  COPY_FIELD(block,ir_,"R");
  COPY_FIELD(block,iv_,"V");
  COPY_FIELD(block,ix_,"X");

  /// refresh Q with callback to p_solver_bicgstab_loop_7 to handle re-entry

  loop_8(block);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_8(EnzoBlock* block) throw() {

  /// access field container on this block

  Field field = block->data()->field();

  if (index_precon_ >= 0) {

    enzo_float* Y = (enzo_float*) field.values(iy_);

    for (int i=0; i<m_; i++) Y[i] = 0.0;

    /// Access the preconditioner for this solver, if any

    Solver * precon = cello::solver(index_precon_);

    /// Apply the preconditioner, then return to
    /// p_solver_bicgstab_loop_8()
    
    precon->set_sync_id (10);
    precon->set_callback(CkIndex_EnzoBlock::p_solver_bicgstab_loop_8());

    /// LINE 10: Y = M \ Q
    
#ifdef TRACE_SOLVER_BCG    
    CkPrintf ("%s %s:%d TRACE_SOLVER_BCG calling preconditioner\n",
	      block->name().c_str(),__FILE__,__LINE__);
#endif    

#ifdef DEBUG_WRITE
    char buffer[40];
    sprintf (buffer,"Q.%s",block->name().c_str());
    FILE * fp = fopen(buffer,"w");
    enzo_float * Q = (enzo_float*) field.values(iq_);
    for (int i=0; i<m_; i++) fprintf (fp,"%20.15g\n",Q[i]);
    fclose(fp);
#endif    

    precon->set_field_x(iy_);
    precon->set_field_b(iq_);

    precon->apply(A_,block);
    
  } else {

    enzo_float * Y = (enzo_float*) field.values(iy_);
    enzo_float * Q = (enzo_float*) field.values(iq_);

    /// LINE 10: Y = M \ Q  [ M = I ]

    for (int i=0; i<m_; i++)  Y[i] = Q[i];

    loop_85(block);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_8() {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_85(this);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_85 (EnzoBlock * block) throw() {
  
  // Refresh field faces then call p_solver_bicgstab_loop_85()

  const int min_face_rank = cello::rank() - 1;
  Refresh refresh (4,min_face_rank,neighbor_type_(), sync_type_(),
		   enzo_sync_id_solver_bicgstab_loop_85);
  

  refresh.set_active(is_finest_(block));

  refresh.add_field (ix_);
  refresh.add_field (ir_);
  refresh.add_field (ir0_);
  refresh.add_field (ip_);
  refresh.add_field (iy_);
  refresh.add_field (iv_);
  refresh.add_field (iq_);
  refresh.add_field (iu_);
  
  block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_bicgstab_loop_9(),&refresh);

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_9() {

  performance_start_(perf_compute,__FILE__,__LINE__);
  
  static_cast<EnzoSolverBiCgStab*> (solver())->loop_10(this);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
  
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_10(EnzoBlock* block) throw() {
  
  /// access field container on this block

  Field field = block->data()->field();

  COPY_FIELD(block,iq_,"Q2_bcg");
  COPY_FIELD(block,iy_,"Y2_bcg");

  if (is_finest_(block)) {

    /// LINE 11:     U = A * Y
    
    A_->matvec(iu_, iy_, block);     /// apply matrix to local block

  }

  COPY_FIELD(block,iu_,"U");

  std::vector<long double> reduce;
  reduce.resize(5+TEMP_INC);
  reduce.clear();
  if (TEMP_INC) reduce[0] = 5;
  
  if (is_finest_(block)) {
    
    enzo_float* U  = (enzo_float*) field.values(iu_);
    enzo_float* Q  = (enzo_float*) field.values(iq_);
    
    /// omega_n = DOT(U, Q)
    /// omega_d = DOT(U, U)
  
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[0+TEMP_INC] += U[i]*Q[i];
	  reduce[1+TEMP_INC] += U[i]*U[i];
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
	    reduce[2+TEMP_INC] += Y[i];
	    reduce[3+TEMP_INC] += U[i];
	    reduce[4+TEMP_INC] += Q[i];
	  }
	}
      }
    }
  }
  
  /// compute sums over Blocks and continue with r_solver_bicgstab_loop_11()

  std::vector<int> is_array;
  if (solve_type_ == solve_tree) {
    is_array.resize(5);
    is_array[0] = is_omega_n_;
    is_array[1] = is_omega_d_;
    is_array[2] = is_ys_;
    is_array[3] = is_us_;
    is_array[4] = is_qs_;
  }

  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_solver_bicgstab_loop_11(NULL), 
     block->proxy_array());
  
  inner_product_(block,5,&reduce[0],is_array,callback);
    
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_11(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_12(this,msg);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_12(EnzoBlock* block,
				 CkReductionMsg * msg) throw() {

  long double* data = (long double*) msg->getData();
  
  double omega_n = data[0+TEMP_INC];
  double omega_d = data[1+TEMP_INC];
  double ys      = data[2+TEMP_INC];
  double us      = data[3+TEMP_INC];
  double qs      = data[4+TEMP_INC];

  delete msg;

  /// verify legal floating-point values for preceding reduction results

  cello::check(omega_n,"BCG_omega_n_",__FILE__,__LINE__);
  cello::check(omega_d,"BCG_omega_d_",__FILE__,__LINE__);

  /// access field container on this block

  Field field = block->data()->field();

  /// for singular problems, update omega_d and project Y and U into R(A)

  if (A_->is_singular()) {

    omega_n -= us*qs/c_;
    omega_d -= us*us/c_;

    if (is_finest_(block)) {
      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* U = (enzo_float*) field.values(iu_);
      enzo_float y_shift = ys / c_;
      enzo_float u_shift = us / c_;
      for (int i=0; i<m_; i++) {
	Y[i] -= y_shift;
	U[i] -= u_shift;
      }
    }
  }

  /// avoid division by 0.0
  
  if (omega_d == 0.0)  omega_d = 1.0;
  
  /// LINE 12:     omega = (U*Q) / (U*U)
  
  omega_ = omega_n / omega_d;

  /// check for breakdown in BiCgStab
  
  if ( omega_ == 0.0 ) {
    WARNING ("EnzoSolverBiCgStab::loop12()",
	     "Solver error: omega_ == 0");
    this->end(block, return_error);
  }

  /// update vectors on leaf blocks
  
  if (is_finest_(block)) {

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
  TRACE_SCALAR(block,"beta_d_",beta_d_);
  
  /// rr_     = DOT(R, R)
  /// beta_n_ = DOT(R, R0)
  
  std::vector<long double> reduce;
  reduce.resize(2+TEMP_INC);
  reduce.clear();
  if (TEMP_INC) reduce[0] = 2;
  
  if (is_finest_(block)) {
    
    enzo_float* R  = (enzo_float*) field.values(ir_);
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[0+TEMP_INC] += R[i]*R[i];
	  reduce[1+TEMP_INC] += R[i]*R0[i];
	}
      }
    }
  }

  /// sum over blocks and continue with r_solver_bicgstab_loop_13()
  
  std::vector<int> is_array;
  if (solve_type_ == solve_tree) {
    is_array.resize(2);
    is_array[0] = is_rr_;
    is_array[1] = is_beta_n_;
  }

  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_solver_bicgstab_loop_13(NULL), 
     block->proxy_array());

  inner_product_(block,2,&reduce[0],is_array,callback);
    
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_13(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_14(this,msg);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_14(EnzoBlock* block,
				 CkReductionMsg * msg) throw() {

  long double* data = (long double*) msg->getData();
  
  rr_     = data[0+TEMP_INC];
  beta_n_ = data[1+TEMP_INC];
  
  TRACE_SCALAR(block,"rr_",rr_);
  TRACE_SCALAR(block,"beta_n_",beta_n_);
  delete msg;

  /// verify legal floating-point values for preceding reduction results
  
  cello::check(rr_,    "BCG_rr_",    __FILE__,__LINE__);
  cello::check(beta_n_,"BCG_beta_n_",__FILE__,__LINE__);

  /// access field container on this block

  Field field = block->data()->field();

  /// check for breakdown in BiCgStab
  
  if (beta_n_ == 0.0) {
    WARNING ("EnzoSolverBiCgStab::loop14()",
	     "Solver error: beta_n_ == 0");
    this->end(block, return_error);
  }
  

  /// LINE 15: beta = (R*R0) / beta_n * (alpha/omega)
  
  enzo_float beta = (beta_n_/beta_d_)*(alpha_/omega_);

  /// update direction vector
  
  if (is_finest_(block)) {

    enzo_float* P = (enzo_float*) field.values(ip_);
    enzo_float* R = (enzo_float*) field.values(ir_);
    enzo_float* V = (enzo_float*) field.values(iv_);

    /// LINE 16:     P = R + beta * (P - omega * V)

    for (int i=0; i<m_; i++) {
      P[i] = R[i] + beta*(P[i] - omega_*V[i]);
    }
  }

  /// contribute to global iteration counter and continue with
  /// r_solver_bicgstab_loop_15()
  
  int iter = iter_ + 1;

  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_solver_bicgstab_loop_15(NULL), 
     block->proxy_array());

  block->contribute(sizeof(int), &iter, 
		    CkReduction::max_int, callback);

}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_15(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_0b(this,msg);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::end (EnzoBlock* block, int retval) throw () {


  deallocate_temporary_(block);
  
  Solver::end_(block);
}

//======================================================================

void EnzoSolverBiCgStab::inner_product_
(EnzoBlock * block, int n, long double * reduce,
 const std::vector<int> & is_array,
 CkCallback callback)
{
#ifdef NEW_CONTRIBUTE
  if (solve_type_ == solve_tree) {
    dot_compute_tree_(block,n,reduce,is_array,callback);
  } else {
    block->contribute((n+1)*sizeof(long double), reduce, 
		      sum_long_double_n_type, callback);
  }
#else
  const int length = n*sizeof(long double);
  switch (n) {
  case 1: 
    block->contribute(length, reduce, sum_long_double_type, callback);
    break;
  case 2:
    block->contribute(length, reduce, sum_long_double_2_type, callback);
    break;
  case 3:
    block->contribute(length, reduce, sum_long_double_3_type, callback);
    break;
  case 4:
    block->contribute(length, reduce, sum_long_double_4_type, callback);
    break;
  case 5:
    block->contribute(length, reduce, sum_long_double_5_type, callback);
    break;
  default:
    ERROR1("EnzoSolverBiCgStab::inner_product_()",
	   "Unsupported n = %d", n);
  }
#endif
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_compute_tree_(EnzoBlock * block,
					   int n,
					   long double * dot_local,
					   const std::vector<int> & is_array,
					   CkCallback callback)
{
  const int level = block->level();

  dot_clear_(block,n,is_array);
  if (level < 0) {
    for (int i=0; i<n; i++) dot_local[i] = 0.0;
    dot_copy_(block,n, dot_local,is_array);
    dot_done_(block,callback);
  } else if (is_finest_(block)) {
    if (level > 0) {
      dot_send_parent_(block,n,dot_local,is_array,callback);
    } else {
      dot_copy_(block,n, dot_local, is_array);
      dot_done_(block,callback);
    }
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_send_parent_(EnzoBlock * block,
					  int n,
					  long double * dot_block,
					  const std::vector<int> & is_array,
					  CkCallback callback)
{
  ASSERT1("EnzoSolverBiCgStab::dot_send_parent()",
	  "level %d must be > 0",
	  block->level(), (block->level() > 0));

  Index index_parent = block->index().index_parent(min_level_);
  
  enzo::block_array()[index_parent].p_dot_recv_parent(n,dot_block,
						      is_array,callback);

}

//----------------------------------------------------------------------

void EnzoBlock::p_dot_recv_parent(int n, long double * dot_block,
				  std::vector<int> is_array,
				  CkCallback callback)
{
  
  //      scalar_dot_partial += dot_child;
  //      if (sync_parent.next()) {
  //         dot_block = scalar_dot_partial;
  //	 scalar_dot_partial = 0.0;
  //         if (level > 0) {
  //            dot_send_parent(dot_block)
  //         } else {
  //	    assert (level == 0);
  //            dot_send_children(dot_block)
  //	    done();
  //         }
  //      }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_send_children_(EnzoBlock * block,
					    int n,
					    long double * dot_local,
					    const std::vector<int> & is_array,
					    CkCallback callback)
{
      // assert(not leaf);
      // for (index_child) {
      //    p_dot_recv_children[index_child](dot_tree)
      // }
}

//----------------------------------------------------------------------

void EnzoBlock::p_dot_recv_children(int n, long double * dot_block,
				    std::vector<int> is_array,
				    CkCallback callback)
{
  //      if (leaf) {
  //         dot_final = dot_tree;
  //      } else {
  //         dot_send_children(dot_tree)
  //      }
  //      dot_done();
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_copy_
(EnzoBlock * block,int n, long double * data, const std::vector<int> & is_array)
{
  Scalar<long double> scalar =
    block->data()->scalar_long_double();
  
  for (int i=0; i<n; i++) {
    *(scalar.value(is_array[i])) = data[i];
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_clear_
(EnzoBlock * block,int n, const std::vector<int> & is_array)
{
  Scalar<long double> scalar =
    block->data()->scalar_long_double();
  
  for (int i=0; i<n; i++) {
    *(scalar.value(is_array[i])) = 0.0;
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_done_(EnzoBlock * block,
				   CkCallback callback)
{
  // call callback_ with block
}

//----------------------------------------------------------------------
