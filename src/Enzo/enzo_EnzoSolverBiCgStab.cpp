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

// #define DEBUG_CALLBACK

// #define TRACE_DOT
// #define TRACE_DOT_CYCLE 0

// #define TRACE_BCG

// #define DEBUG_FIELD
// #define DEBUG_REDUCE
// #define DEBUG_GHOST
// #define DEBUG_SCALAR

// #define DEBUG_WRITE
// #define DEBUG_READ

#define S(index) scalar_(block,is_##index##_)

#ifdef DEBUG_SCALAR
#   define TRACE_SCALAR(BLOCK,NAME,SCALAR)				\
  CkPrintf ("%s:%d %s TRACE_SCALAR %s = %Lg\n",				\
	    __FILE__,__LINE__,BLOCK->name().c_str(),NAME,SCALAR);
#else
#   define TRACE_SCALAR(BLOCK,NAME,SCALAR) /* ... */
#endif

#ifdef DEBUG_GHOST
#   define TRACE_GHOST(BLOCK,ID,NAME)					\
  {									\
    Field field = BLOCK->data()->field();				\
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
	      " %Lg %Lg %Lg\n"						\
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
    long double sum_a=0.0,sum_abs=0.0;					\
    for (int iz=gz_; iz<mz_-gz_; iz++) {				\
      for (int iy=gy_; iy<my_-gy_; iy++) {				\
	for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	  int i=ix+mx_*(iy+my_*iz);					\
	  sum_a+=X[i];							\
	  sum_abs+=std::abs(X[i]);					\
	}								\
      }									\
    }									\
    CkPrintf ("%s:%d %s %s COPY_FIELD %d %s shift %Lg %Lg\n"		\
	      ,__FILE__,__LINE__,BLOCK->name().c_str(),name().c_str(),ID,COPY,sum_a, sum_abs); \
    TRACE_GHOST(BLOCK,ID,COPY);						\
  }
#else
#   define COPY_FIELD(BLOCK,ID,COPY) /* ... */
#endif

#ifdef TRACE_BCG
#  undef TRACE_BCG
#  define TRACE_BCG(BLOCK,SOLVER,msg)					\
  if (BLOCK->cycle() >= TRACE_DOT_CYCLE) {				\
    CkPrintf ("%d %s %s:%d TRACE_BCG %s %s level %d\n",			\
	      CkMyPe(), (BLOCK!=NULL) ? BLOCK->name().c_str():"no block", \
	      __FILE__,__LINE__,					\
	      (SOLVER!=NULL) ? SOLVER->name().c_str():"no solver",msg,	\
	      (BLOCK!=NULL) ? BLOCK->level() : -99);			\
    fflush(stdout);							\
  }
#else
#  define TRACE_BCG(BLOCK,SOLVER,msg) /* ... */
#endif
  
#ifdef TRACE_DOT
#   undef TRACE_DOT
#  define TRACE_DOT(BLOCK,msg,ifunction)				\
  if (BLOCK->cycle() >= TRACE_DOT_CYCLE) {				\
    CkPrintf ("%d %s %s:%d TRACE_DOT %s %d\n",				\
	      CkMyPe(), (BLOCK!=NULL) ? BLOCK->name().c_str():"no block", \
	      __FILE__,__LINE__,					\
	      msg, ifunction);						\
    fflush(stdout);							\
  }
#else
#  define TRACE_DOT(BLOCK,msg,ifunction) /* ... */
#endif


//----------------------------------------------------------------------

EnzoSolverBiCgStab::EnzoSolverBiCgStab
(std::string name,
 std::string field_x, std::string field_b,
 int monitor_iter, int restart_cycle,
 int solve_type,
 int min_level, int max_level,
 int iter_max, double res_tol,
 int index_precon,
 int coarse_level
 ) 
  : Solver(name,
	   field_x,
	   field_b,
	   monitor_iter,
	   restart_cycle,
	   solve_type,
	   min_level,
	   max_level),
    res_tol_(res_tol),
    A_(NULL),
    index_precon_(index_precon),
    iter_max_(iter_max), 
    ir_(0), ir0_(0), ip_(0), 
    iy_(0), iv_(0), iq_(0), iu_(0),
    m_(0), mx_(0), my_(0), mz_(0),
    gx_(0), gy_(0), gz_(0),
    coarse_level_(coarse_level)
{

  //  if (solve_type == solve_tree) {
  ScalarDescr * scalar_descr_quad = cello::scalar_descr_long_double();

  // skip index==0 for checking index validity  
  int is_skip   =  scalar_descr_quad->new_value("solver_bicgstab_skip");
  
  is_alpha_ =  scalar_descr_quad->new_value("solver_bicgstab_alpha");
  is_beta_n_ = scalar_descr_quad->new_value("solver_bicgstab_beta_n");
  is_beta_d_ = scalar_descr_quad->new_value("solver_bicgstab_beta_d");
  is_omega_ =  scalar_descr_quad->new_value("solver_bicgstab_omega");
  is_omega_d_ =scalar_descr_quad->new_value("solver_bicgstab_omega_d");
  is_omega_n_ =scalar_descr_quad->new_value("solver_bicgstab_omega_n");
  is_err_ =    scalar_descr_quad->new_value("solver_bicgstab_err");
  is_err0_ =   scalar_descr_quad->new_value("solver_bicgstab_err0");
  is_err_min_ =scalar_descr_quad->new_value("solver_bicgstab_err_min");
  is_err_max_ =scalar_descr_quad->new_value("solver_bicgstab_err_max");
  is_rr_ =     scalar_descr_quad->new_value("solver_bicgstab_rr");
  is_r0s_ =    scalar_descr_quad->new_value("solver_bicgstab_r0s");
  is_c_ =      scalar_descr_quad->new_value("solver_bicgstab_c");
  is_bs_ =     scalar_descr_quad->new_value("solver_bicgstab_bs");
  is_xs_ =     scalar_descr_quad->new_value("solver_bicgstab_xs");
  is_bnorm_ =  scalar_descr_quad->new_value("solver_bicgstab_bnorm");
  is_rho0_  =  scalar_descr_quad->new_value("solver_bicgstab_rho0");
  is_vr0_ =    scalar_descr_quad->new_value("solver_bicgstab_vr0");
  is_ys_ =     scalar_descr_quad->new_value("solver_bicgstab_ys");
  is_vs_ =     scalar_descr_quad->new_value("solver_bicgstab_vs");
  is_us_ =     scalar_descr_quad->new_value("solver_bicgstab_us");
  is_qs_ =     scalar_descr_quad->new_value("solver_bicgstab_qs");

  if (solve_type == solve_tree) {
   
    function_.push_back(&EnzoSolverBiCgStab::start_2); // inner_product 0
    function_.push_back(&EnzoSolverBiCgStab::loop_0a); // inner_product 1
    function_.push_back(&EnzoSolverBiCgStab::loop_6);  // inner_product 2
    function_.push_back(&EnzoSolverBiCgStab::loop_12); // inner_product 3
    function_.push_back(&EnzoSolverBiCgStab::loop_14); // inner_product 4
    
    ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
    is_dot_sync_ = scalar_descr_sync->new_value("solver_bicgstab_dot_sync");

  }
  
  ScalarDescr * scalar_descr_int = cello::scalar_descr_int();
  is_iter_ = scalar_descr_int->new_value("solver_bicgstab_iter");

  FieldDescr * field_descr = cello::field_descr();

  ir_ = field_descr->insert_temporary();
  ir0_ = field_descr->insert_temporary();
  ip_ = field_descr->insert_temporary();
  iy_ = field_descr->insert_temporary();
  iv_ = field_descr->insert_temporary();
  iq_ = field_descr->insert_temporary();
  iu_ = field_descr->insert_temporary();
  
  /// Initialize default Refresh (called before entry to compute())

  // upper limit on A_->ghost_depth(), which isn't known yet
  const int ghost_depth = 4; 
  const int min_face_rank = cello::rank() - 1;

  const int ir = add_refresh
    (ghost_depth,min_face_rank,
     neighbor_type_(), sync_type_(),
     enzo_sync_id_solver_bicgstab);

  if (solve_type_ == solve_tree)
    refresh(ir) -> set_root_level (coarse_level_);

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

  TRACE_BCG(block,this,"apply");
  
  Solver::begin_(block);
  
  EnzoBlock* enzo_block = enzo::block(block);

  if (solve_type_ == solve_tree) {
    s_dot_sync_(enzo_block) = cello::num_children();
  }
  
  A_ = A;

  Field field = block->data()->field();

  allocate_temporary_(block);

  /// cast input argument to the EnzoBlock associated with this char

  /// access the field infromation on this block
  
  field.dimensions (0, &mx_, &my_, &mz_);
  field.ghost_depth(0, &gx_, &gy_, &gz_);

  m_ = mx_*my_*mz_;

  compute_ (enzo_block);
}

//======================================================================

void EnzoSolverBiCgStab::compute_(EnzoBlock* block) throw() {

  TRACE_BCG(block,this,"compute");
  
  /// initialize BiCgStab iteration counter
  (s_iter_(block)) = 0;

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

  COPY_FIELD(block,ib_,"B0_bcg");
  for (int i=0; i<m_; i++) {
    X[i] = R[i] = R0[i] = P[i] = 0.0;
    Y[i] = V[i] = Q[i] =  U[i] = 0.0;
  }

  if (is_finest_(block)) {

    const bool reuse_x = reuse_solution_ (block->cycle());
    
    if ( reuse_x ) {

      enzo_float* X_copy  = (enzo_float*) field.values("X_copy");

      for (int i=0; i<m_; i++) X[i] = X_copy[i];

      A_->residual (ir_, ib_, ix_, block);
      
    }
  }

  /// for singular Poisson problems, N(A) is not empty, so project B
  /// into R(A)

  if (is_singular_()) {

    std::vector<long double> reduce;
    reduce.resize(3+1);
    reduce.clear();
    reduce[0] = 3;

    long double count = 0.0;

    if (is_finest_(block)) {

      enzo_float* B = (enzo_float*) field.values(ib_);
      enzo_float* X = (enzo_float*) field.values(ix_);
 
      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    count++;
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[2] += B[i];
	    reduce[3] += X[i];
	  }
	}
      }
      reduce[1] = count;
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

    CkCallback callback
      (CkIndex_EnzoBlock::r_solver_bicgstab_start_1(NULL),
       block->proxy_array());

#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::r_solver_bicgstab_start_1(NULL));
#endif    

    TRACE_DOT(block,"start",0);
    inner_product_(block,3,&reduce[0],is_array,callback,0); // start_2

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


  TRACE_BCG(block,this,"start_2");
  
  if (solve_type_ != solve_tree && msg != NULL) {
    long double* data = (long double*) msg->getData();
    ASSERT1("EnzoSolverBiCgStab::start_2",
	    "Expecting (data[0] = %d) == 3",
	    data[0],(data[0] == 3));
    S(c)  = data[1];
    S(bs) = data[2];
    S(xs) = data[3];
  }

  delete msg;
  
  /// access field container on this block

  Field field = block->data()->field();

  /// update B and initialize temporary vectors (on leaf blocks only)

  std::vector<long double> reduce;
  reduce.resize(3+1);
  reduce.clear();
  reduce[0] = 3;

  if (is_finest_(block)) {

    /// access relevant fields
    enzo_float* B  = (enzo_float*) field.values(ib_);
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    enzo_float* P  = (enzo_float*) field.values(ip_);
    enzo_float* R  = (enzo_float*) field.values(ir_);
    enzo_float* X  = (enzo_float*) field.values(ix_);

    /// for singular problems, project B into R(A)

    if (is_singular_()) {

      long double bs = S(bs);
      long double xs = S(xs);
      long double c  = S(c);

      enzo_float b_shift = bs / c;
      enzo_float x_shift = xs / c;

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
	  reduce[1] += R[i]*R0[i];
	  reduce[2] += B[i]*B[i];
	  reduce[3] += R[i];
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

  CkCallback callback (CkIndex_EnzoBlock::r_solver_bicgstab_start_3(NULL), 
		       block->proxy_array());
#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::r_solver_bicgstab_start_3(NULL));
#endif    


  TRACE_DOT(block,"start",1);
  inner_product_(block,3,&reduce[0],is_array,callback,1); // loop_0a

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

  TRACE_BCG(block,this,"loop_0a");
  
  if (solve_type_ != solve_tree && msg != NULL) {
    long double* data = (long double*) msg->getData();
    ASSERT1("EnzoSolverBiCgStab::loop_0a",
	    "Expecting (data[0] = %d) == 3",
	    data[0],(data[0] == 3));
    S(beta_n) = data[1];
    S(bnorm)  = data[2];
    S(r0s)    = data[3];
  }

  delete msg;

  TRACE_SCALAR(block,"beta_n",S(beta_n));
  TRACE_SCALAR(block,"bnorm_",S(bnorm));
  TRACE_SCALAR(block,"r0s_",S(r0s));

  loop_0(block);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_0b(EnzoBlock* block,
				 CkReductionMsg * msg) throw() {

  delete msg;

  TRACE_BCG(block,this,"loop_0b");

  (s_iter_(block))++;
  loop_0(block);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_0(EnzoBlock* block) throw() {

  /// verify legal floating-point value for preceding reduction result

  TRACE_BCG(block,this,"loop_0");

  if (is_finest_(block)) cello::check(S(beta_n),"BCG_beta_n",__FILE__,__LINE__);

  /// initialize/update current error, store error statistics
  
  const int cycle = block->cycle();
  const bool reuse_x = reuse_solution_ (cycle);

  const int iter = (s_iter_(block));
  if (iter == 0) {
    const long double s_r0s = S(r0s);
    const long double s_c =   S(c);
    if (is_singular_()) {
      enzo_float shift = s_r0s/ s_c;
	
      if (is_finest_(block)) {
	Field field = block->data()->field();
	enzo_float* R0 = (enzo_float*) field.values(ir0_);
	enzo_float* P  = (enzo_float*) field.values(ip_);
	enzo_float* R  = (enzo_float*) field.values(ir_);
	for (int i=0; i<m_; i++) {
	  R0[i] -= shift;
	  R[i]  -= shift;
	  P[i]  -= shift;
	}
	S(beta_n) -= s_r0s*s_r0s/s_c;
	TRACE_SCALAR(block,"beta_n",S(beta_n));
      }
    }
    S(rho0) = sqrt(S(bnorm)); // ||B||
    
    TRACE_SCALAR(block,"rho0_",S(rho0));
    if (S(rho0) == 0.0)  S(rho0) = 1.0;
    S(err) =
      sqrt(S(beta_n)) / S(rho0);
    S(err0)    = S(err);
    S(err_min) = S(err);
    S(err_max) = S(err);
  } else {
    S(err) =
      sqrt(S(rr)) / S(rho0);
    S(err_min) =
      std::min(S(err), S(err_min));
    S(err_max) =
      std::max(S(err), S(err_max));
  }
  TRACE_SCALAR(block,"err_",S(err));


  const bool is_converged = (S(err) < res_tol_);
  const bool is_diverged  = (iter >= iter_max_);

  /// monitor output solution progress (iteration, residual, etc)

  int a3[3];
  block->index().array(a3,a3+1,a3+2);
  
  const bool l_output =
    ( (a3[0]==0 && a3[1]==0 && a3[2]==0) &&
      (block->level()==coarse_level_) &&
      ( (iter == 0) ||
	(is_converged || is_diverged) ||
	(monitor_iter_ && (iter % monitor_iter_) == 0 )) );

  if (l_output) {
    if (reuse_x) {
      monitor_output_(block,iter,sqrt(S(bnorm)),
		      S(err_min)*sqrt(S(bnorm)),
		      S(err)    *sqrt(S(bnorm)),
		      S(err_max)*sqrt(S(bnorm)));
    } else {
      monitor_output_(block,iter,
		      S(err0),
		      S(err_min),
		      S(err),
		      S(err_max));
    }
  }

  /// Write final status if done
  if (block->index().is_root() && (is_converged || is_diverged) ) {
    CkPrintf ("%s DEBUG_SOLVER bicgstab "
	      "final iter = %d rr = %Lg  rho0 = %Lg  rr/rho0 = %Lg\n",
	      block->name().c_str(),
	      iter,S(rr),S(rho0),sqrt(S(rr))/ S(rho0));
  }

  if (is_converged) {

    /// Save copy of X if it will be used as initial guess next cycle

    const bool reuse_next_x = reuse_solution_ (cycle + 1);

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

  TRACE_BCG(block,this,"loop_2");

  Field field = block->data()->field();

  if (index_precon_ >= 0) {

    enzo_float * Y = (enzo_float*) field.values(iy_);

    for (int i=0; i<m_; i++) Y[i] = 0.0;

    /// Access the preconditioner for this solver, if any

    Solver * precon = cello::solver(index_precon_);

    /// Apply the preconditioner, then return to
    /// p_solver_bicgstab_loop_2()

    precon->set_sync_id (enzo_sync_id_solver_bicgstab_precon_1);
    precon->set_callback(CkIndex_EnzoBlock::p_solver_bicgstab_loop_2());

#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::p_solver_bicgstab_loop_2());
#endif    

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
    
    

    precon->set_field_x(iy_);
    precon->set_field_b(ip_);
    precon->apply(A_,block);
    
  } else { // no preconditioner

    enzo_float * Y = (enzo_float*) field.values(iy_);
    enzo_float * P = (enzo_float*) field.values(ip_);
    
    /// LINE 04: Y = M \ P  [ M = I ]
    for (int i=0; i<m_; i++) Y[i] = P[i];

    loop_25(block);

  }
  COPY_FIELD(block,ip_,"P0_bcg");
  COPY_FIELD(block,iy_,"Y0_bcg");
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_bicgstab_loop_2() {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_25(this);
  
  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::loop_25 (EnzoBlock * block) throw() {
  
  TRACE_BCG(block,this,"loop_25");

  // Refresh field faces then call p_solver_bicgstab_loop_25()
  
  const int ghost_depth = A_->ghost_depth();
  const int min_face_rank = cello::rank() - 1;

  Refresh refresh
    (ghost_depth,min_face_rank,
     neighbor_type_(), sync_type_(),
     enzo_sync_id_solver_bicgstab_loop_25);

  if (solve_type_ == solve_tree)
    refresh.set_root_level (coarse_level_);

  refresh.set_active(is_finest_(block));

  refresh.add_field (ix_);

  refresh.add_field (ir_);
  refresh.add_field (ir0_);
  refresh.add_field (ip_);
  refresh.add_field (iy_);
  refresh.add_field (iv_);
  refresh.add_field (iq_);
  refresh.add_field (iu_);

#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::p_solver_bicgstab_loop_3());
#endif    
    
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

  TRACE_BCG(block,this,"loop_4");
  
  /// access field container on this block

  Field field = block->data()->field();

  /// V = MATVEC(A,Y)

  COPY_FIELD(block,iy_,"Y1_bcg");
  COPY_FIELD(block,ip_,"P1_bcg");
  
  if (is_finest_(block)) {

    /// LINE 05: V = A * Y
    
    A_->matvec(iv_, iy_, block);

  }

  COPY_FIELD(block,iv_,"V1_bcg");
  /// compute local contributions to vr0_ = DOT(V, R0)
  
  std::vector<long double> reduce;
  reduce.resize(3+1);
  reduce.clear();
  reduce[0] = 3;
  
  if (is_finest_(block)) {
    
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    enzo_float* V  = (enzo_float*) field.values(iv_);

    /// LINE 07 [part]  vr0_ = V*R0
    
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[1] += V[i]*R0[i];
	}
      }
    }
    
    /// for singular Poisson problems need all vectors in R(A), so
    /// project both Y and V into R(A)

    if (is_singular_()) {

      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* V = (enzo_float*) field.values(iv_);

      /// ys_ = sum (Y[i])
      /// vs_ = sum (V[i])
      
      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[2] += Y[i];
	    reduce[3] += V[i];
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

#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::r_solver_bicgstab_loop_5(NULL));
#endif    
  TRACE_DOT(block,"start",2);
  inner_product_(block,3,&reduce[0],is_array,callback,2); // loop_6
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

  TRACE_BCG(block,this,"loop_6");

  if (solve_type_ != solve_tree && msg != NULL) {
    long double* data = (long double*) msg->getData();
    ASSERT1("EnzoSolverBiCgStab::loop_6",
	    "Expecting (data[0] = %d) == 3",
	    data[0],(data[0] == 3));
    S(vr0) = data[1];
    S(ys)  = data[2];
    S(vs)  = data[3];
  }

  delete msg;

  const long double vr0 = S(vr0);
  const long double ys =  S(ys);
  const long double vs =  S(vs);
  
  TRACE_SCALAR(block,"vr0",vr0);
  TRACE_SCALAR(block,"ys",ys);
  TRACE_SCALAR(block,"vs",vs);

  COPY_FIELD(block,iv_,"V1_bcg");

  Field field = block->data()->field();

  if (is_finest_(block)) {

    /// for singular problems, project Y and V into R(A)

    if (is_singular_()) {
    
      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* V = (enzo_float*) field.values(iv_);

      enzo_float y_shift = ys / S(c);
      enzo_float v_shift = vs / S(c);

      for (int i=0; i<m_; i++) {
	Y[i] -= y_shift;
	V[i] -= v_shift;
      }
      COPY_FIELD(block,iy_,"Y_shift");
      COPY_FIELD(block,iv_,"V_shift");

    }

    /// compute alpha factor in BiCgStab algorithm (all blocks)

    /// LINE 07:     alpha = beta_n / (V*R0)

    S(alpha) = S(beta_n) / vr0;

    TRACE_SCALAR(block,"alpha",S(alpha));

    /// update vectors (on leaf blocks only)

    /// access relevant fields
    enzo_float* Q = (enzo_float*) field.values(iq_);
    enzo_float* R = (enzo_float*) field.values(ir_);
    enzo_float* V = (enzo_float*) field.values(iv_);
    enzo_float* X = (enzo_float*) field.values(ix_);
    enzo_float* Y = (enzo_float*) field.values(iy_);

    /// LINE 08: Q = R - alpha * V
    /// LINE 09: X = X + alpha * Y
    
    COPY_FIELD(block,ir_,"R1_bcg");
    enzo_float alpha = S(alpha);
    
    for (int i=0; i<m_; i++) {
      Q[i] = R[i] - alpha*V[i];
      X[i] = X[i] + alpha*Y[i];
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

  TRACE_BCG(block,this,"loop_8");

  /// access field container on this block

  Field field = block->data()->field();

  if (index_precon_ >= 0) {

    enzo_float* Y = (enzo_float*) field.values(iy_);

    for (int i=0; i<m_; i++) Y[i] = 0.0;

    /// Access the preconditioner for this solver, if any

    Solver * precon = cello::solver(index_precon_);

    /// Apply the preconditioner, then return to
    /// p_solver_bicgstab_loop_8()
    
    precon->set_sync_id (enzo_sync_id_solver_bicgstab_precon_2);
    precon->set_callback(CkIndex_EnzoBlock::p_solver_bicgstab_loop_8());
#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::p_solver_bicgstab_loop_8());
#endif    

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
    for (int i=0; i<m_; i++) fprintf (fp,"%g\n",Q[i]);
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
  
  TRACE_BCG(block,this,"loop_85");

  // Refresh field faces then call p_solver_bicgstab_loop_85()

  const int ghost_depth = A_->ghost_depth();
  const int min_face_rank = cello::rank() - 1;

  Refresh refresh
    (ghost_depth,min_face_rank,
     neighbor_type_(), sync_type_(),
     enzo_sync_id_solver_bicgstab_loop_85);
  
  if (solve_type_ == solve_tree)
    refresh.set_root_level (coarse_level_);

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

  TRACE_BCG(block,this,"loop_10");

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
  reduce.resize(5+1);
  reduce.clear();
  reduce[0] = 5;
  
  if (is_finest_(block)) {
    
    enzo_float* U  = (enzo_float*) field.values(iu_);
    enzo_float* Q  = (enzo_float*) field.values(iq_);
    
    /// omega_n = DOT(U, Q)
    /// omega_d = DOT(U, U)
  
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[1] += U[i]*Q[i];
	  reduce[2] += U[i]*U[i];
	}
      }
    }
  
    /// for singular Poisson problems, project both Y and U into R(A)

    if (is_singular_()) {

      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* U = (enzo_float*) field.values(iu_);

      /// ys_ = SUM(Y)
      /// us_ = SUM(U)

      for (int iz=gz_; iz<mz_-gz_; iz++) {
	for (int iy=gy_; iy<my_-gy_; iy++) {
	  for (int ix=gx_; ix<mx_-gx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[3] += Y[i];
	    reduce[4] += U[i];
	    reduce[5] += Q[i];
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

#ifdef DEBUG_REDUCE  
  CkPrintf ("DEBUG_REDUCE %s %s:%d %Lg %Lg %Lg %Lg %Lg\n",
	    block->name().c_str(),__FILE__,__LINE__,
	    reduce[0],reduce[1],reduce[2],reduce[3],reduce[4]);
#endif
  
  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_solver_bicgstab_loop_11(NULL), 
     block->proxy_array());
  
#ifdef DEBUG_CALLBACK    
  CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	    __FILE__,__LINE__,CkIndex_EnzoBlock::r_solver_bicgstab_loop_11(NULL));
#endif    

  TRACE_DOT(block,"start",3);
  inner_product_(block,5,&reduce[0],is_array,callback,3); // loop_12
    
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

  TRACE_BCG(block,this,"loop_12");

  if (solve_type_ != solve_tree && msg != NULL) {
    long double* data = (long double*) msg->getData();
    ASSERT1("EnzoSolverBiCgStab::loop_12",
	    "Expecting (data[0] = %d) == 5",
	    data[0],(data[0] == 5));
    S(omega_n) = data[1];
    S(omega_d) = data[2];
    S(ys)      = data[3];
    S(us)      = data[4];
    S(qs)      = data[5];
  }

  delete msg;

  const long double ys = S(ys);
  const long double us = S(us);
  const long double qs = S(qs);

  /// verify legal floating-point values for preceding reduction results

  if (is_finest_(block)) {
    cello::check(S(omega_n),"BCG_omega_n_",__FILE__,__LINE__);
    cello::check(S(omega_d),"BCG_omega_d_",__FILE__,__LINE__);
  }

  /// access field container on this block

  Field field = block->data()->field();

  /// for singular problems, update omega_d and project Y and U into R(A)

  if (is_singular_()) {

    S(omega_n) -= us*qs/ S(c);
    S(omega_d) -= us*us/ S(c);

    if (is_finest_(block)) {
      enzo_float* Y = (enzo_float*) field.values(iy_);
      enzo_float* U = (enzo_float*) field.values(iu_);
      enzo_float y_shift = ys / S(c);
      enzo_float u_shift = us / S(c);
      for (int i=0; i<m_; i++) {
	Y[i] -= y_shift;
	U[i] -= u_shift;
      }
    }
  }
  
  /// avoid division by 0.0
  
  if (S(omega_d) == 0.0)  S(omega_d) = 1.0;
  
  /// LINE 12:     omega = (U*Q) / (U*U)
  
  S(omega) = S(omega_n) / S(omega_d);

  /// check for breakdown in BiCgStab
  
  if ( S(omega_n) == 0.0 ) {
    WARNING1 ("EnzoSolverBiCgStab::loop12()",
	      "Solver error: %s omega_n == 0",
	      block->name().c_str());
    this->end(block, return_error);
  }
  if ( S(omega_d) == 0.0 ) {
    WARNING1 ("EnzoSolverBiCgStab::loop12()",
	     "Solver error: %s omega_d1 == 0",
	      block->name().c_str());
    this->end(block, return_error);
  }
  if ( S(omega) == 0.0 ) {
    WARNING1 ("EnzoSolverBiCgStab::loop12()",
	     "Solver error: %s omega_ == 0",
	      block->name().c_str());
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
      X[i] = X[i] + S(omega)*Y[i];
      R[i] = Q[i] - S(omega)*U[i];
    }
  }

  
  /// Update previous beta value (beta_d_) to current value (beta_n_)
  
  S(beta_d) = S(beta_n);
  
  /// rr_     = DOT(R, R)
  /// beta_n = DOT(R, R0)
  
  std::vector<long double> reduce;
  reduce.resize(2+1);
  reduce.clear();
  reduce[0] = 2;
  
  if (is_finest_(block)) {
    
    enzo_float* R  = (enzo_float*) field.values(ir_);
    enzo_float* R0 = (enzo_float*) field.values(ir0_);
    
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  reduce[1] += R[i]*R[i];
	  reduce[2] += R[i]*R0[i];
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

#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::r_solver_bicgstab_loop_13(NULL));
#endif    

  TRACE_DOT(block,"start",4);
  inner_product_(block,2,&reduce[0],is_array,callback,4); // loop_14
    
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

  TRACE_BCG(block,this,"loop_14");

  if (solve_type_ != solve_tree && msg != NULL) {
    long double* data = (long double*) msg->getData();
    ASSERT1("EnzoSolverBiCgStab::loop_14",
	    "Expecting (data[0] = %d) == 2",
	    data[0],(data[0] == 2));
    S(rr)     = data[1];
    S(beta_n) = data[2];
  }
  
  delete msg;

  TRACE_SCALAR(block,"rr_",S(rr));
  TRACE_SCALAR(block,"beta_n_",S(beta_n));

  /// verify legal floating-point values for preceding reduction results

  if (is_finest_(block)) {
    cello::check(S(rr),    "BCG_rr_",   __FILE__,__LINE__);
    cello::check(S(beta_n),"BCG_beta_n",__FILE__,__LINE__);
  }

  /// access field container on this block

  Field field = block->data()->field();

  /// check for breakdown in BiCgStab
  
  if (S(beta_n) == 0.0) {
    WARNING1 ("EnzoSolverBiCgStab::loop14()",
	     "Solver error: %s beta_n == 0",
	      block->name().c_str());
    this->end(block, return_error);
  }
  

  /// LINE 15: beta = (R*R0) / beta_n * (alpha/omega)

  enzo_float beta = (S(beta_n)/S(beta_d))*(S(alpha)/ S(omega));

  /// update direction vector
  
  if (is_finest_(block)) {

    enzo_float* P = (enzo_float*) field.values(ip_);
    enzo_float* R = (enzo_float*) field.values(ir_);
    enzo_float* V = (enzo_float*) field.values(iv_);

    /// LINE 16:     P = R + beta * (P - omega * V)

    for (int i=0; i<m_; i++) {
      P[i] = R[i] + beta*(P[i] - S(omega)*V[i]);
    }
  }

  /// contribute to global iteration counter and continue with
  /// r_solver_bicgstab_loop_15()
  
  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_solver_bicgstab_loop_15(NULL), 
     block->proxy_array());

#ifdef DEBUG_CALLBACK    
    CkPrintf ("DEBUG_CALLBACK %s:%d %d\n",
	      __FILE__,__LINE__,CkIndex_EnzoBlock::r_solver_bicgstab_loop_15(NULL));
#endif    

  loop_0b(block,NULL);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_bicgstab_loop_15(CkReductionMsg* msg) {

  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverBiCgStab*> (solver())->loop_0b(this,msg);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::end (EnzoBlock* block, int retval) throw () {

  TRACE_BCG(block,this,"end");

  deallocate_temporary_(block);
  
  Solver::end_(block);
}

//======================================================================

void EnzoSolverBiCgStab::inner_product_
(EnzoBlock * block, int n, long double * reduce,
 const std::vector<int> & is_array,
 CkCallback callback,
 int i_function)
{
  if (solve_type_ == solve_tree) {
    TRACE_BCG(block,this,"inner_product_A");
    dot_compute_tree_(block,n,reduce+1,is_array,i_function);
  } else {
    TRACE_BCG(block,this,"inner_product_B");
    block->contribute((n+1)*sizeof(long double), reduce, 
		      sum_long_double_n_type, callback);
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_compute_tree_(EnzoBlock * block,
					   int n,
					   long double * dot_local,
					   const std::vector<int> & is_array,
					   int i_function)
{
  TRACE_DOT(block,"dot_compute_tree",i_function);
  dot_clear_(block,n,is_array);

  const int level = block->level();
  if (level < coarse_level_) {
    dot_done_(block,i_function,__FILE__,__LINE__);
  } else if (is_finest_(block)) {
    if (level > coarse_level_) {
      dot_send_parent_(block,n,dot_local,is_array,i_function);
    } else {
      dot_save_(block,n, dot_local, is_array);
      dot_done_(block,i_function,__FILE__,__LINE__);
    }
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_send_parent_(EnzoBlock * block,
					  int n,
					  long double * dot_block,
					  const std::vector<int> & is_array,
					  int i_function)
{
  TRACE_DOT(block,"dot_send_parent",i_function);
  ASSERT2("EnzoSolverBiCgStab::dot_send_parent()",
	  "level %d must be > coarse_level = %d",
	  block->level(), coarse_level_,
	  (block->level() > coarse_level_));

  Index index_parent = block->index().index_parent(min_level_);

  enzo::block_array()[index_parent].p_dot_recv_parent(n,dot_block,
						      is_array,i_function);

}

//----------------------------------------------------------------------

void EnzoBlock::p_dot_recv_parent(int n, long double * dot_block,
				  std::vector<int> is_array,
				  int i_function)
{
  auto solver = static_cast<EnzoSolverBiCgStab*> (this->solver());

  solver->dot_recv_parent(this,n,dot_block,is_array,i_function);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_recv_parent(EnzoBlock * block,
					 int n,
					 long double * dot_block,
					 const std::vector<int> & is_array,
					 int i_function)
{
  
  TRACE_DOT(block,"dot_recv_parent",i_function);
  dot_increment_(block,n,is_array,dot_block);
  
  Sync & sync = s_dot_sync_(block);
  if (sync.next()) {
    dot_load_(block,n, dot_block, is_array);
    if (block->level() > coarse_level_) {
      //      dot_clear_(block,n,is_array);
      dot_send_parent_(block,n,dot_block,is_array,i_function);
    } else {
      dot_send_children_(block,n,dot_block,is_array,i_function);
      dot_done_(block,i_function,__FILE__,__LINE__);
    }
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_send_children_(EnzoBlock * block,
					    int n,
					    long double * dot_local,
					    const std::vector<int> & is_array,
					    int i_function)
{
  TRACE_DOT(block,"dot_send_children",i_function);
  ItChild it_child(cello::rank());
  int ic3[3];
  while (it_child.next(ic3)) {

    Index index_child = block->index().index_child(ic3,min_level_);

    enzo::block_array()[index_child].p_dot_recv_children
      (n,dot_local,is_array,i_function);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_dot_recv_children(int n, long double * dot_block,
				    std::vector<int> is_array,
				    int i_function)
{
  auto solver = static_cast<EnzoSolverBiCgStab*> (this->solver());
  solver->dot_recv_children(this,n,dot_block,is_array,i_function);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_recv_children(EnzoBlock * block,
					   int n,
					   long double * dot_local,
					   const std::vector<int> & is_array,
					   int i_function)
{
  TRACE_DOT(block,"dot_recv_children",i_function);
  dot_save_(block,n, dot_local, is_array);
  if (!is_finest_(block)) {
    dot_send_children_(block,n,dot_local,is_array,i_function);
  }
  dot_done_(block,i_function,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_save_
(EnzoBlock * block,int n, long double * data, const std::vector<int> & is_array)
{
  TRACE_DOT(block,"dot_save",-1);
  Scalar<long double> scalar =
    block->data()->scalar_long_double();
  
  for (int i=0; i<n; i++) {
    *(scalar.value(is_array[i])) = data[i];
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_load_
(EnzoBlock * block,int n, long double * data, const std::vector<int> & is_array)
{
  TRACE_DOT(block,"dot_load",-1);
  Scalar<long double> scalar =
    block->data()->scalar_long_double();
  
  for (int i=0; i<n; i++) {
    data[i] = *(scalar.value(is_array[i]));
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

void EnzoSolverBiCgStab::dot_clear_
(EnzoBlock * block,int n, long double * array)
{
  for (int i=0; i<n; i++) array[i] = 0.0;
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_increment_
(EnzoBlock * block,
 int n,
 const std::vector<int> & is_array,
 long double * dot_block)
{
  TRACE_DOT(block,"dot_increment",-1);
  Scalar<long double> scalar =
    block->data()->scalar_long_double();
  
  for (int i=0; i<n; i++) {
    *(scalar.value(is_array[i])) += dot_block[i];
  }
}

//----------------------------------------------------------------------

void EnzoSolverBiCgStab::dot_done_(EnzoBlock * block,
				   int i_function,
				   const char * file, int line)
{
  TRACE_DOT(block,"dot_done",i_function);
  (this->*function_[i_function])(block,NULL);
}

//----------------------------------------------------------------------
