// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityBiCGStab.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     2014-10-23 16:19:06
/// @brief    Implements the EnzoMethodGravityBiCGStab class
//
// The following is Matlab code for right-preconditioned BiCGStab, 
// optimized for both memory efficiency and for reduced 'reduction' 
// synchronization:
//
// function [x, err, iter, flag] = bcgs_opt(A, b, M, iter_max, res_tol, singular)
// % Usage: [x, err, iter, flag] = bcgs_opt(A, b, M, iter_max, res_tol, singular)
// %
// % Solves the linear system Ax=b using the BiConjugate Gradient
// % Stabilized Method with right preconditioning.  We assume an
// % initial guess of x=0.
// %
// % input   A        REAL matrix
// %         b        REAL right hand side vector
// %         M        REAL preconditioner matrix
// %         iter_max INTEGER maximum number of iterations
// %         res_tol  REAL error tolerance
// %         singular INTEGER flag denoting singular matrix
// %
// % output  x        REAL solution vector
// %         err      REAL error norm
// %         iter     INTEGER number of iterations performed
// %         flag     INTEGER: 0 = solution found to tolerance
// %                           1 = no convergence given iter_max
// %                          -1 = breakdown: rho = 0
// %                          -2 = breakdown: omega = 0
//
//   % initialize outputs
//   iter = 0;
//   flag = 0;
//
//   % for singular (Poisson) systems, project RHS to range space
//   if (singular)
//      bs = sum(b);                    % global reduction
//      bc = length(b);                 % global reduction
//      b = b - bs/bc;
//   end
//   
//   % compute initial residual, initialize temporary vectors
//   r = b;
//   r0 = r;
//   p = r;                             % Note: r, r0 and p are in R(A)
//   
//   % inner-product at initialization: <b,b>
//   rho0 = dot(b,b);                   % global reduction
//   beta_n = rho0;
// 
//   % store ||b|| for relative error checks
//   rho0 = sqrt(rho0);
//   if (rho0 == 0.0), rho0 = 1.0; end
//
//   % check initial relative residual to see if work is finished
//   err = sqrt(beta_n)/rho0;
//   if (err < res_tol), return; end
//
//   % begin iteration
//   for iter = 1:iter_max
//
//      % first half of step
//      y = M\p;                        % ? communication
//      v = A*y;                        % pt-to-pt communication
//      vr0 = dot(v,r0);                % global reduction
//      if (singular)
//         ys = sum(y);                 % global reduction
//         vs = sum(v);                 % global reduction
//      end
//      if (singular)                % project y and v into R(A)
//         y = y - ys/bc;            % note: fine to do after vr0 since <r0,e>=0
//         v = v - vs/bc;
//      end
//      alpha = beta_n / vr0;
//      q = r - alpha*v;
//      x = x + alpha*y;             % Note: q and x are in R(A)
//
//      % stabilization portion of step
//      y = M\q;                        % ? communication
//      u = A*y;                        % pt-to-pt communication
//      omega_d = dot(u,u);             % global reduction
//      omega_n = dot(u,q);             % global reduction
//      if (singular)
//         ys = sum(y);                 % global reduction
//         us = sum(u);                 % global reduction
//      end
//      if (singular)                % project y and u into R(A)
//         y = y - ys/bc;            % note: fine to do after omega_n since <q,e>=0
//         u = u - us/bc;
//         omega_d=omega_d-us^2/bc;  % fix due to reduction before projection
//      end
//      if (omega_d == 0),  omega_d = 1; end
//      omega = omega_n / omega_d;
//      if (omega == 0), flag = -2; return; end
//      x = x + omega*y;
//      r = q - omega*u;             % Note: r and x are in R(A)
//      if (beta_n == 0), flag = -1; return; end
//      err = sqrt(rr)/rho0;
//      if (err < res_tol), return; end
//
//      % compute new direction vector
//      beta = (beta_n/beta_d)*(alpha/omega);
//      p = r + beta*(p - omega*v);  % Note: p is in R(A)
//
//   end      
//
//   % non-convergent
//   flag = 1;
// 
// % END bcgs.m
//
//     V -- (temporary) matrix-vector product vector
//     Q -- (temporary) temporary vector
//     U -- (temporary) temporary vector
//   
//   Communicated scalars (results of inner-products):
//     rho0 -- floating-point, initial residual
//     beta_d -- floating-point
//     beta_n -- floating-point
//     vr0 -- floating-point
//     omega_d -- floating-point
//     omega_n -- floating-point
//     rr -- floating-point
//     bs -- floating-point, shift factor numerator
//     bc -- floating-point, shift factor denominator
//
//   Local scalars:
//     iter -- integer
//     beta -- floating-point
//     err -- floating-point
//     err_min -- floating-point
//     err_max -- floating-point
//     alpha -- floating-point
//     omega -- floating-point
//     iter_max -- integer (input)
//     res_tol -- floating-point (input)
//
// ======================================================================
//
// BiCGStab partitioned along parallel communication / synchronization steps
//
// --------------------
// compute_()
// --------------------
//
//    return_ = return_unknown;
// 
//    B = <right-hand side>
//    X = <initial solution X0> => initialize to zero
//    
//    iter = 0
//
//    if (is_singular_) {
//       bs_ = SUM(B)   ==> r_gravity_bicgstab_start_1
//       bc_ = COUNT(B) ==> r_gravity_bicgstab_start_1
//    } else {
//       call gravity_bicgstab_start_2
//    }  
//
// --------------------
// r_gravity_bicgstab_start_1()
// --------------------
//
//    receive bs_ and bc_
//
//    call start_2
//
// --------------------
// start_2()
// --------------------
//
//    if (is_singular_)
//      B = B - bs_/bc_;
//
//    R = B
//    R0 = R
//    P = R
//    beta_n_ = DOT(R, R) ==> r_gravity_bicgstab_start_3
//
// --------------------
// r_gravity_bicgstab_start_3()
// --------------------
//
//    receive beta_n_
// 
//    call loop_0
//
// ==================================================
//
// --------------------
// loop_0()
// --------------------
//
//    if (iter_ == 0) {
//       rho0_ = sqrt(beta_n_)
//       if (rho0_ == 0.0)  rho0_ = 1.0
//       err_ = sqrt(beta_n_) / rho0_
//       err_min_ = err_
//       err_max_ = err_
//    } else {
//       err_ = sqrt(rr_) / rho0_;
//       err_min_ = min(err_, err_min_)
//       err_max_ = max(err_, err_max_)
//    }
//
//    output solution progress (iteration, residual, etc)
//
//    if (err_ < res_tol_)
//       ==> end(return_converged)
//
//    if (iter_ >= iter_max_) {
//       ==> end(return_error_max_iter_reached);
//    }
// 
//    refresh(P) ==> p_gravity_bicgstab_loop_1
//
// --------------------
// p_gravity_bicgstab_loop_1()
// --------------------
//
//    call loop_2
//
// --------------------
// loop_2()
// --------------------
//
//    Y = SOLVE(M,P)
//
//    refresh(Y) ==> p_gravity_bicgstab_loop_3
//
// --------------------
// p_gravity_bicgstab_loop_3()
// --------------------
//  
//    call loop_4
// 
// --------------------
// loop_4()
// --------------------
//
//    V = MATVEC(A,Y)
//
//    vr0_ = DOT(V, R0) ==> r_gravity_bicgstab_loop_5
//
//    if (is_singular_) {
//       ys_ = SUM(Y) ==> r_gravity_bicgstab_loop_5
//       vs_ = SUM(V) ==> r_gravity_bicgstab_loop_5
//    }  
//
// --------------------
// r_gravity_bicgstab_loop_5()
// --------------------
//  
//    receive vr0_, ys_ and vs_
//
//    call loop_6
// 
// --------------------
// loop_6()
// --------------------
//
//    if (is_singular_) {
//      Y = Y - ys_/bc_;
//      V = V - vs_/bc_;
//    }
//
//    alpha_ = beta_n_ / ( vr0_ );
//    Q = R - alpha_*V
//    X = X + alpha_*Y
//
//    refresh(Q) ==> p_gravity_bicgstab_loop_8
//
// --------------------
// p_gravity_bicgstab_loop_7()
// --------------------
//  
//    call loop_8
//
// --------------------
// loop_8()
// --------------------
//  
//    Y = SOLVE(M,Q)
//
//    refresh(Y) ==> p_gravity_bicgstab_loop_9
// 
// --------------------
// p_gravity_bicgstab_loop_9()
// --------------------
//  
//    call loop_10
// 
// --------------------
// loop_10()
// --------------------
//
//    U = MATVEC(A,Y)
//
//    omega_d_ = DOT(U, U) ==> r_gravity_bicgstab_loop_11
//    omega_n_ = DOT(U, Q) ==> r_gravity_bicgstab_loop_11
//
//    if (is_singular_) {
//       ys_ = SUM(Y) ==> r_gravity_bicgstab_loop_11
//       us_ = SUM(U) ==> r_gravity_bicgstab_loop_11
//    }  
//
// --------------------
// r_gravity_bicgstab_loop_11()
// --------------------
//  
//    receive ys_, us_, omega_d_ and omega_n_
//
//    call loop_12
// 
// --------------------
// loop_12()
// --------------------
//
//    if (is_singular_) {
//      Y = Y - ys_/bc_;
//      U = U - us_/bc_;
//
//      omega_d_ = omega_d - us_*us_/bc;
//    }
//
//    if (omega_d_ == 0.0)  omega_d_ = 1.0
//    if ( omega_n_ == 0.0 ) 
//       end(return_error_omega_eq_0)
//
//    omega_ = omega_n_ / omega_d_
//    X = X + omega_*Y
//    R = Q - omega_*U
//
//    beta_d_ = beta_n_
//
//    rr_     = DOT(R, R)  ==> r_gravity_bicgstab_loop_13
//    beta_n_ = DOT(R, R0) ==> r_gravity_bicgstab_loop_13
//
// --------------------
// r_gravity_bicgstab_loop_13()
// --------------------
//
//    receive rr_ and beta_n_
//
//    call loop_14
//
// --------------------
// loop_14()
// --------------------
//
//    if ( beta_n_ == 0.0 )
//       end(return_error_beta_n_eq_0)
//
//    beta_ = (beta_n_/beta_d_)*(alpha_/omega_)
//    P = R + beta_*( P - omega_*V )
//
//    iter = iter_ + 1
//    (contribute iter to iter_)
//
//    ==> r_gravity_bicgstab_loop_15()
//
// --------------------
// r_gravity_bicgstab_loop_15()
// --------------------
//
//    receive iter_
//
//    call loop_0
//
//
// ==================================================
//
// --------------------
// end(return_)
// --------------------
//
//    if (retval == return_converged) {
//       potential = X
//       compute acceleration field
//       ==> gravity_bicgstab_exit()
//    } else {
//       ERROR (retval)
//    }
//
// --------------------------------------------------


#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

//----------------------------------------------------------------------

EnzoMethodGravityBiCGStab::EnzoMethodGravityBiCGStab
(FieldDescr* field_descr, int rank,
 double grav_const, int iter_max, double res_tol, 
 int monitor_iter, bool is_singular, bool diag_precon) 
  : Method(), 
    A_(new EnzoMatrixLaplace),
    M_(NULL),
    is_singular_(is_singular),
    rank_(rank),
    grav_const_(grav_const),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    monitor_iter_(monitor_iter),
    rho0_(0), err_(0), err_min_(0), err_max_(0),
    idensity_(0),  ipotential_(0),
    ib_(0), ix_(0), ir_(0), ir0_(0), ip_(0), 
    iy_(0), iv_(0), iq_(0), iu_(0),
    nx_(0), ny_(0), nz_(0),
    mx_(0), my_(0), mz_(0),
    gx_(0), gy_(0), gz_(0),
    iter_(0),
    beta_d_(0), beta_n_(0), beta_(0), 
    omega_d_(0), omega_n_(0), omega_(0), 
    vr0_(0), rr_(0), alpha_(0),
    bs_(0.0),bc_(0.0),
    ys_(0.0),vs_(0.0),us_(0.0),
    id_refresh_P_(-1),
    id_refresh_Q_(-1),
    id_refresh_Y_(-1)
{

  /// set preconditioning operator
  M_ = (diag_precon) ? (Matrix*)(new EnzoMatrixDiagonal) 
                     : (Matrix*)(new EnzoMatrixIdentity);

  /// access problem-defining fields for eventual RHS and solution
  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");

  /// access existing fields for temporary vectors (must be declared in parameter file)
  ib_ = field_descr->field_id(name() + "_B");
  ix_ = field_descr->field_id(name() + "_X");
  ir_ = field_descr->field_id(name() + "_R");
  ir0_ = field_descr->field_id(name() + "_R0");
  ip_ = field_descr->field_id(name() + "_P");
  iy_ = field_descr->field_id(name() + "_Y");
  iv_ = field_descr->field_id(name() + "_V");
  iq_ = field_descr->field_id(name() + "_Q");
  iu_ = field_descr->field_id(name() + "_U");

  // /// create temporary fields for use within solver
  // int precision = field_descr->precision(idensity_);
  // ib_ = field_descr->insert_temporary(name() + "_B");
  // field_descr->set_precision(ib_, precision);
  // ix_ = field_descr->insert_temporary(name() + "_X");
  // field_descr->set_precision(ix_, precision);
  // ir_ = field_descr->insert_temporary(name() + "_R");
  // field_descr->set_precision(ir_, precision);
  // ir0_ = field_descr->insert_temporary(name() + "_R0");
  // field_descr->set_precision(ir0_, precision);
  // ip_ = field_descr->insert_temporary(name() + "_P");
  // field_descr->set_precision(ip_, precision);
  // iy_ = field_descr->insert_temporary(name() + "_Y");
  // field_descr->set_precision(iy_, precision);
  // iv_ = field_descr->insert_temporary(name() + "_V");
  // field_descr->set_precision(iv_, precision);
  // iq_ = field_descr->insert_temporary(name() + "_Q");
  // field_descr->set_precision(iq_, precision);
  // iu_ = field_descr->insert_temporary(name() + "_U");
  // field_descr->set_precision(iu_, precision);


  /// Initialize default Refresh (called before entry to compute())
  const int num_fields = field_descr->field_count();  /// DR: NECESSARY?
  const int ir = add_refresh(1, rank-1, neighbor_leaf, sync_barrier);
  refresh(ir)->add_field(idensity_);

  /// Initialize specific vector Refreshes
  id_refresh_P_ = add_refresh(1, rank-1, neighbor_leaf, sync_barrier);
  refresh(id_refresh_P_)->add_field(ip_);
  id_refresh_Q_ = add_refresh(1, rank-1, neighbor_leaf, sync_barrier);
  refresh(id_refresh_Q_)->add_field(iq_);
  id_refresh_Y_ = add_refresh(1, rank-1, neighbor_leaf, sync_barrier);
  refresh(id_refresh_Y_)->add_field(iy_);

}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::compute(Block* block) throw() {

  /// cast input argument to the EnzoBlock associated with this char
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (block);

  /// access the field infromation on this block
  Field field = block->data()->field();
  field.size(&nx_, &ny_, &nz_);
  field.dimensions(idensity_, &mx_, &my_, &mz_);
  field.ghost_depth(idensity_, &gx_, &gy_, &gz_);

  // /// allocate temporary vector data for use witin solve
  // field.allocate_temporary(ib_);
  // field.allocate_temporary(ix_);
  // field.allocate_temporary(ir_);
  // field.allocate_temporary(ir0_);
  // field.allocate_temporary(ip_);
  // field.allocate_temporary(iy_);
  // field.allocate_temporary(iv_);
  // field.allocate_temporary(iq_);
  // field.allocate_temporary(iu_);

  /// call templated internal compute_ routine
  int precision = field.precision(idensity_);
  if      (precision == precision_single)    compute_<float>      (enzo_block);
  else if (precision == precision_double)    compute_<double>     (enzo_block);
  else if (precision == precision_quadruple) compute_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
}

//======================================================================

/// eventually replace this with a more flexible reduction type
/// (this one requires an array with 4 entries but handles long doubles)
extern CkReduction::reducerType r_method_gravity_bicgstab_type;

template<class T> void EnzoMethodGravityBiCGStab::compute_(EnzoBlock* enzo_block) throw() {

  /// initialize BiCGStab iteration counter
  iter_ = 0;

  /// access field container on this block
  Data* data  = enzo_block->data();
  Field field = data->field();

  /// construct RHS B, initialize initial solution X to zero (only on leaf blocks)
  if (enzo_block->is_leaf()) {

    /// access relevant fields
    T* density = (T*) field.values(idensity_);
    T* B       = (T*) field.values(ib_);
    T* X       = (T*) field.values(ix_);

    /// set X = 0 [Q: necessary?  couldn't we reuse the solution from the previous solve?]
    /// set B = -h^2 * 4 * PI * G * density
    for (int iz=0; iz<mz_; iz++) 
      for (int iy=0; iy<my_; iy++) 
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  X[i] = 0.0;
	  B[i] = - 4.0 * (cello::pi) * grav_const_ * density[i];
	}
  }

  /// for singular Poisson problems, N(A) is not empty, so project B into R(A)
  if (is_singular_) {

    /// set bs_ = SUM(B)   ==> r_gravity_bicgstab_start_1
    /// set bc_ = COUNT(B) ==> r_gravity_bicgstab_start_1
    long double reduce[4] = {0.0, 0.0, 0.0, 0.0};
    if (enzo_block->is_leaf()) {
      T* B = (T*) field.values(ib_);
      reduce[0] = sum_(B);
      reduce[1] = count_();
    }

    /// initiate callback for r_gravity_bicgstab_start_1 and contribute to sum and count
    CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_start_1<T>(NULL), 
			enzo_block->proxy_array());
    enzo_block->contribute(4*sizeof(long double), &reduce, 
			   r_method_gravity_bicgstab_type, callback);

  } else {

    /// nonsingular system, just call start_2 directly
    this->start_2<T>(enzo_block);

  }
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_start_1(CkReductionMsg* msg) {

  /// EnzoBlock accumulates global contributions to SUM(B) and COUNT(B)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  long double* data = (long double*) msg->getData();
  method->set_bs( data[0] );
  method->set_bc( data[1] );
  delete msg;

  /// call start_2 to continue
  method->start_2<T>(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::start_2(EnzoBlock* enzo_block) throw() {

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();

  /// update B and initialize temporary vectors (on leaf blocks only)
  if (enzo_block->is_leaf()) {

    /// access relevant fields
    T* B  = (T*) field.values(ib_);
    T* R0 = (T*) field.values(ir0_);
    T* P  = (T*) field.values(ip_);
    T* R  = (T*) field.values(ir_);

    /// for singular problems, project B into R(A)
    if (is_singular_) {
      T shift = -bs_ / bc_;
      for (int iz=0; iz<mz_; iz++) 
	for (int iy=0; iy<my_; iy++) 
	  for (int ix=0; ix<mx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    B[i] += shift;
	  }
    }

    /// initialize R = R0 = P = B
    for (int iz=0; iz<mz_; iz++) 
      for (int iy=0; iy<my_; iy++) 
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  R[i] = R0[i] = P[i] = B[i];
	}

  }

  /// Compute local contributions to beta_n_ = DOT(R, R)
  long double reduce[4] = {0.0, 0.0, 0.0, 0.0};
  if (enzo_block->is_leaf()) {
    T* R = (T*) field.values(ir_);
    reduce[0] = dot_(R,R);
  }
  
  /// initiate callback for r_gravity_bicgstab_start_3 and contribute to dot-product
  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_start_3<T>(NULL), 
		      enzo_block->proxy_array());
  enzo_block->contribute(4*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_start_3(CkReductionMsg* msg) {

  /// EnzoBlock accumulates global contributions to DOT(R, R)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  long double* data = (long double*) msg->getData();
  method->set_beta_n( data[0] );
  delete msg;

  /// call loop_0 to begin solver loop
  method->loop_0<T>(this);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_0(EnzoBlock* enzo_block) throw() {

  /// verify legal floating-point value for preceding reduction result
  cello::check(beta_n_,"beta_n_",__FILE__,__LINE__);

  /// initialize/update current error, store error statistics
  if (iter_ == 0) {
    rho0_ = sqrt(beta_n_);
    if (rho0_ == 0.0)  rho0_ = 1.0;
    err_ = sqrt(beta_n_) / rho0_;
    err_min_ = err_;
    err_max_ = err_;
  } else {
    err_ = sqrt(rr_) / rho0_;
    err_min_ = std::min(err_, err_min_);
    err_max_ = std::max(err_, err_max_);
  }

  /// output solution progress (iteration, residual, etc)
  if (enzo_block->index().is_root()) {
    Monitor* monitor = enzo_block->simulation()->monitor();
    if (iter_ == 0)  
      monitor->print("Enzo", "BiCGStab iter %04d  rho0 %.16g",
		     iter_,(double)(rho0_));
    if (monitor_iter_ && (iter_ % monitor_iter_) == 0 ) 
      monitor_output_(enzo_block);
  }

  /// check for convergence
  if (err_ < res_tol_) {
    this->end<T>(enzo_block, return_converged);

  /// otherwise check for failure
  } else if (iter_ >= iter_max_)  {
    this->end<T>(enzo_block, return_error_max_iter_reached);

  /// otherwise refresh P with callback to p_gravity_bicgstab_loop_1 for re-entry
  } else {

    this->refresh(id_refresh_P_)->set_active(enzo_block->is_leaf());
    enzo_block->refresh_enter(CkIndex_EnzoBlock::p_gravity_bicgstab_loop_1(),
			      this->refresh(id_refresh_P_));

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_gravity_bicgstab_loop_1() {
  /// re-entry into loop_2, using template with appropriate precision
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  int precision = field.precision(field.field_id("density"));  // assumes field precisions equal
  if      (precision == precision_single)    
    method->loop_2<float>(enzo_block);
  else if (precision == precision_double)    
    method->loop_2<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->loop_2<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_2(EnzoBlock* enzo_block) throw() {

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;

  /// Y = SOLVE(M, P)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);  /// set mesh increments (used in M_->matvec)
    M_->matvec(iy_, ip_, enzo_block);     /// apply preconditioner to local block
  }

  // refresh Y with callback to p_gravity_bicgstab_loop_3 to handle re-entry
  this->refresh(id_refresh_Y_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_enter(CkIndex_EnzoBlock::p_gravity_bicgstab_loop_3(),
			    this->refresh(id_refresh_Y_));

}

//----------------------------------------------------------------------

void EnzoBlock::p_gravity_bicgstab_loop_3() {
  /// re-entry into loop_4, using template with appropriate precision
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  if      (precision == precision_single)    
    method->loop_4<float>(enzo_block);
  else if (precision == precision_double)    
    method->loop_4<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->loop_4<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_4(EnzoBlock* enzo_block) throw() {

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;

  /// V = MATVEC(A,Y)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);  /// set mesh increments (used in A_->matvec)
    A_->matvec(iv_, iy_, enzo_block);     /// apply matrix to local block
  }

  /// compute local contributions to vr0_ = DOT(V, R0)
  long double reduce[4] = {0.0, 0.0, 0.0, 0.0};
  if (enzo_block->is_leaf()) {
    T* R0 = (T*) field.values(ir0_);
    T* V  = (T*) field.values(iv_);
    reduce[0] = dot_(V,R0);
  }

  /// for singular Poisson problems need all vectors in R(A), so project both Y and V into R(A)
  if (is_singular_) {

    /// set ys_ = SUM(Y)
    /// set vs_ = SUM(V)
    if (enzo_block->is_leaf()) {
      T* Y = (T*) field.values(iy_);
      T* V = (T*) field.values(iv_);
      reduce[1] = sum_(Y);
      reduce[2] = sum_(V);
    }

  }

  /// initiate callback to r_gravity_bicgstab_loop_5 and contribute to global sums
  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_5<T>(NULL), 
		      enzo_block->proxy_array());
  enzo_block->contribute(4*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);


}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_5(CkReductionMsg* msg) {

  /// EnzoBlock accumulates global contributions to SUM(Y) and SUM(V)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  long double* data = (long double*) msg->getData();
  method->set_vr0( data[0] );
  method->set_ys(  data[1] );
  method->set_vs(  data[2] );
  delete msg;

  /// call loop_6 to continue
  method->loop_6<T>(this);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_6(EnzoBlock* enzo_block) throw() {

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();

  /// for singular problems, project Y and V into R(A)
  if (enzo_block->is_leaf() && is_singular_) {
    T* Y = (T*) field.values(iy_);
    T* V = (T*) field.values(iv_);
    T yshift = -ys_ / bc_;
    T vshift = -vs_ / bc_;
    for (int iz=0; iz<mz_; iz++) 
      for (int iy=0; iy<my_; iy++) 
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  Y[i] += yshift;
	  V[i] += vshift;
	}
  }

  /// compute alpha factor in BiCGStab algorithm (all blocks)
  alpha_ = beta_n_ / vr0_;

  /// update vectors (on leaf blocks only)
  if (enzo_block->is_leaf()) {

    /// access relevant fields
    T* Q = (T*) field.values(iq_);
    T* R = (T*) field.values(ir_);
    T* V = (T*) field.values(iv_);
    T* X = (T*) field.values(ix_);
    T* Y = (T*) field.values(iy_);

    /// update: Q = -alpha_*V + R
    zaxpy_(Q, -alpha_, V, R);

    /// update: X = alpha_*Y + X
    zaxpy_(X, alpha_, Y, X);

  }

  /// refresh Q with callback to p_gravity_bicgstab_loop_7 to handle re-entry
  this->refresh(id_refresh_Q_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_enter(CkIndex_EnzoBlock::p_gravity_bicgstab_loop_7(),
			    this->refresh(id_refresh_Q_));

}

//----------------------------------------------------------------------

void EnzoBlock::p_gravity_bicgstab_loop_7() {
  /// re-entry into loop_8, using template with appropriate precision
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  if      (precision == precision_single)    
    method->loop_8<float>(enzo_block);
  else if (precision == precision_double)    
    method->loop_8<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->loop_8<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_8(EnzoBlock* enzo_block) throw() {

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;
 
  // Y = SOLVE(M, Q)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);  /// set mesh increments (used in A_->matvec)
    M_->matvec(iy_, iq_, enzo_block);     /// apply preconditioner to local block
  }

  // refresh Y with callback to p_gravity_bicgstab_loop_9 to handle re-entry
  this->refresh(id_refresh_Y_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_enter(CkIndex_EnzoBlock::p_gravity_bicgstab_loop_9(),
			    this->refresh(id_refresh_Y_));

}

//----------------------------------------------------------------------

void EnzoBlock::p_gravity_bicgstab_loop_9() {
  /// re-entry into loop_10, using template with appropriate precision
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  if      (precision == precision_single)    
    method->loop_10<float>(enzo_block);
  else if (precision == precision_double)    
    method->loop_10<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->loop_10<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_10(EnzoBlock* enzo_block) throw() {

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;

  // U = MATVEC(A,Y)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);  /// set mesh increments (used in A_->matvec)
    A_->matvec(iu_, iy_, enzo_block);     /// apply matrix to local block
  }

  /// compute local contributions to
  /// omega_d_ = DOT(U, U)
  /// omega_n_ = DOT(U, Q)
  long double reduce[4] = {0.0, 0.0, 0.0, 0.0};
  if (enzo_block->is_leaf()) {
    T* U  = (T*) field.values(iu_);
    T* Q  = (T*) field.values(iq_);
    reduce[0] = dot_(U,U);
    reduce[1] = dot_(U,Q);
  }

  /// for singular Poisson problems, project both Y and U into R(A)
  if (is_singular_ && enzo_block->is_leaf()) {

    /// set ys_ = SUM(Y)
    /// set us_ = SUM(U)
    T* Y = (T*) field.values(iy_);
    T* U = (T*) field.values(iu_);
    reduce[2] = sum_(Y);
    reduce[3] = sum_(U);
  }

  /// initiate callback to r_gravity_bicgstab_loop_11, and contribute to overall dot-products
  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_11<T>(NULL), 
		      enzo_block->proxy_array());
  enzo_block->contribute(4*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
    
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_11(CkReductionMsg* msg) {

  /// EnzoBlock accumulates global contributions to SUM(Y) and SUM(U)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  long double* data = (long double*) msg->getData();
  method->set_omega_d( data[0] );
  method->set_omega_n( data[1] );
  method->set_ys( data[2] );
  method->set_us( data[3] );
  delete msg;

  /// call loop_12 to continue
  method->loop_12<T>(this);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_12(EnzoBlock* enzo_block) throw() {

  /// verify legal floating-point values for preceding reduction results
  cello::check(omega_d_,"omega_d_",__FILE__,__LINE__);
  cello::check(omega_n_,"omega_n_",__FILE__,__LINE__);

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();

  /// for singular problems, update omega_d_ and project Y and U into R(A)
  if (is_singular_) {

    omega_d_ = omega_d_ - us_*us_/bc_;

    if (enzo_block->is_leaf()) {
      T* Y = (T*) field.values(iy_);
      T* U = (T*) field.values(iu_);
      T yshift = -ys_ / bc_;
      T ushift = -us_ / bc_;
      for (int iz=0; iz<mz_; iz++) 
	for (int iy=0; iy<my_; iy++) 
	  for (int ix=0; ix<mx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    Y[i] += yshift;
	    U[i] += ushift;
	  }
    }
  }

  /// fix omega_d_ if necessary (for division)
  if (omega_d_ == 0.0)  omega_d_ = 1.0;
  
  /// compute omega factor in BiCGStab algorithm (all blocks)
  omega_ = omega_n_ / omega_d_;

  /// check for breakdown in BiCGStab
  if ( omega_ == 0.0 )
    this->end<T>(enzo_block, return_error_omega_eq_0);

  /// update vectors (on leaf blocks only)
  if (enzo_block->is_leaf()) {

    /// access relevant fields
    T* X = (T*) field.values(ix_);
    T* Y = (T*) field.values(iy_);
    T* R = (T*) field.values(ir_);
    T* Q = (T*) field.values(iq_);
    T* U = (T*) field.values(iu_);
    
    /// update: X = omega_*Y + X
    zaxpy_(X, omega_, Y, X);

    /// update: R = -omega_*U + Q
    zaxpy_(R, -omega_, U, Q);

  }

  /// Update previous beta value (beta_d_) to current value (beta_n_)
  beta_d_ = beta_n_;

  /// compute local contributions to
  /// rr_     = DOT(R, R)
  /// beta_n_ = DOT(R, R0)
  long double reduce[4] = {0.0, 0.0, 0.0, 0.0};
  if (enzo_block->is_leaf()) {
    T* R  = (T*) field.values(ir_);
    T* R0 = (T*) field.values(ir0_);
    reduce[0] = dot_(R,R);
    reduce[1] = dot_(R,R0);
  }

  /// initiate callback to r_gravity_bicgstab_loop_13, and contribute to overall dot-products
  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_13<T>(NULL), 
		      enzo_block->proxy_array());
  enzo_block->contribute(4*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
    
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_13(CkReductionMsg* msg) {

  // EnzoBlock accumulates global contributions to DOT(R,R) and DOT(R,R0)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  long double* data = (long double*) msg->getData();
  method->set_rr(     data[0] );
  method->set_beta_n( data[1] );
  delete msg;

  // call loop_14 to continue
  method->loop_14<T>(this);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::loop_14(EnzoBlock* enzo_block) throw() {

  /// verify legal floating-point values for preceding reduction results
  cello::check(rr_,    "rr_",    __FILE__,__LINE__);
  cello::check(beta_n_,"beta_n_",__FILE__,__LINE__);

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();

  /// check for breakdown in BiCGStab
  if (beta_n_ == 0.0)
    this->end<T>(enzo_block, return_error_beta_n_eq_0);
  
  /// compute beta factor in BiCGStab algorithm (all blocks)
  beta_ = (beta_n_/beta_d_)*(alpha_/omega_);

  /// update direction vector (on leaf blocks only) -- P = R+beta*(P-omega*V)
  if (enzo_block->is_leaf()) {

    /// access relevant fields
    T* P = (T*) field.values(ip_);
    T* R = (T*) field.values(ir_);
    T* V = (T*) field.values(iv_);

    /// update: P = beta_*P + R
    zaxpy_(P, beta_, P, R);

    /// update: P = -beta_*omega_*V + P
    zaxpy_(P, -beta_*omega_, V, P);

  }

  /// contribute to global iteration counter
  int iter = iter_ + 1;

  /// initiate callback to r_gravity_bicgstab_loop_15, and contribute to overall max
  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_15<T>(NULL), 
		      enzo_block->proxy_array());
  enzo_block->contribute(sizeof(int), &iter, 
			 CkReduction::max_int, callback);

}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_15(CkReductionMsg* msg) {

  // EnzoBlock accumulates global contributions to iter
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  method->set_iter ( ((int*)msg->getData())[0] );
  delete msg;

  // call loop_0 to continue to next iteration
  method->loop_0<T>(this);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::end(EnzoBlock* enzo_block, int retval) throw () {
  /// supposed to do:
  ///    if (retval == return_converged) {
  ///       potential = X
  ///       compute acceleration field
  ///       ==> gravity_bicgstab_exit()
  ///    } else {
  ///       ERROR (retval)
  ///    }
  /// 
  /// actually just does
  ///    potential = X
  ///    compute acceleration field
  ///    enzo_block->compute_done

  /// access field container on this block
  Data* data = enzo_block->data();
  Field field = data->field();

  /// extract the solution and compute derived acceleration fields (leaf blocks only)
  if (enzo_block->is_leaf()) {

    /// access relevant fields
    T* X         = (T*) field.values(ix_);
    T* potential = (T*) field.values(ipotential_);

    /// output solution progress (iteration, residual, etc)
    if (enzo_block->index().is_root()) 
      monitor_output_(enzo_block);

    /// copy potential from solution vector X
    copy_(potential, X, mx_, my_, mz_);

    /// compute acceleration fields from potential
    bool symmetric;
    int order;
    EnzoComputeAcceleration compute_acceleration(field.field_descr(),
						 rank_, symmetric=true,
						 order=2);
    compute_acceleration.compute(enzo_block);
  }

  // /// deallocate temporary vector data
  // field.deallocate_temporary(ib_);
  // field.deallocate_temporary(ix_);
  // field.deallocate_temporary(ir_);
  // field.deallocate_temporary(ir0_);
  // field.deallocate_temporary(ip_);
  // field.deallocate_temporary(iy_);
  // field.deallocate_temporary(iv_);
  // field.deallocate_temporary(iq_);
  // field.deallocate_temporary(iu_);

  /// indicate that method is finished
  enzo_block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::monitor_output_(EnzoBlock* enzo_block) throw() {

  Monitor* monitor = enzo_block->simulation()->monitor();
  monitor->print("Enzo", "BiCGStab iter %04d  err %.16g [%g %g]",
		 iter_,
		 (double)(err_),
		 (double)(err_min_),
		 (double)(err_max_));
}

//======================================================================

template<class T> long double EnzoMethodGravityBiCGStab::dot_(const T* X, const T* Y) const throw() {

  const int i0 = gx_ + mx_*(gy_ + my_*gz_);
  long double value = 0.0;
  for (int iz=0; iz<nz_; iz++) 
    for (int iy=0; iy<ny_; iy++) 
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + (ix + mx_*(iy + my_*iz));
	value += X[i]*Y[i];
      }
  return value;
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::zaxpy_(T* Z, double a, const T* X, const T* Y) const throw() {

  for (int iz=0; iz<mz_; iz++) 
    for (int iy=0; iy<my_; iy++) 
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	Z[i] = a * X[i] + Y[i];
      }
}

//----------------------------------------------------------------------

template<class T> long double EnzoMethodGravityBiCGStab::sum_(const T* X) const throw() {

  const int i0 = gx_ + mx_*(gy_ + my_*gz_);
  long double value = 0.0;
  for (int iz=0; iz<nz_; iz++) 
    for (int iy=0; iy<ny_; iy++) 
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + (ix + mx_*(iy + my_*iz));
	value += X[i];
      }
  return value;
}

//----------------------------------------------------------------------

long double EnzoMethodGravityBiCGStab::count_() const throw() {
  return 1.0*nx_*ny_*nz_;
}

//----------------------------------------------------------------------
