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
//   p = r;
//   
//   % inner-product at initialization: <b,b>
//   rho0 = dot(b,b);                   % global reduction
//   beta_d = rho0;
// 
//   % store ||b|| for relative error checks
//   rho0 = sqrt(rho0);
//   if (rho0 == 0.0), rho0 = 1.0; end
//
//   % check initial relative residual to see if work is finished
//   err = sqrt(beta_d)/rho0;
//   if (err < res_tol), return; end
//
//   % begin iteration
//   for iter = 1:iter_max
//
//      % first half of step
//      y = M\p;                        % ? communication
//      v = A*y;                        % pt-to-pt communication
//      vr0 = dot(v,r0);                % global reduction
//      alpha = beta_d / vr0;
//      q = r - alpha*v;
//      x = x + alpha*y;
//
//      % stabilization portion of step
//      y = M\q;                        % ? communication
//      u = A*y;                        % pt-to-pt communication
//      omega_d = dot(u,u);             % global reduction
//      omega_n = dot(u,q);             % global reduction
//      if (omega_d == 0),  omega_d = 1; end
//      omega = omega_n / omega_d;
//      if (omega == 0), flag = -2; return; end
//      x = x + omega*y;
//      r = q - omega*u;
//
//      % compute remaining dot products for step
//      rr = dot(r,r);                  % global reduction
//      beta_n = dot(r,r0);             % global reduction
//      
//      % check for failure/convergence
//      if (beta_n == 0), flag = -1; return; end
//      err = sqrt(rr)/rho0;
//      if (err < res_tol), return; end
//
//      % compute new direction vector
//      beta = (beta_n/beta_d)*(alpha/omega);
//      beta_d = beta_n;
//      p = r + beta*(p - omega*v);
//
//   end      
//
//   % non-convergent
//   flag = 1;
// 
// % END bcgs.m
//
// ======================================================================
//
// Data requirements:
//
//   Matrices/linear operators:
//     A -- matrix
//     M -- preconditioner
//
//   Vectors:
//     X -- solution vector
//     B -- right-hand side vector
//     R -- (temporary) residual vector
//     R0 -- (temporary) original residual vector
//     P -- (temporary) direction vector
//     Y -- (temporary) preconditioner solution vector
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
//    bs_ = SUM(B)   ==> r_gravity_bicgstab_start_1
//    bc_ = COUNT(B) ==> r_gravity_bicgstab_start_1
//
// --------------------
// r_gravity_bicgstab_start_1()
// --------------------
//
//    receive bs_ and bc_
//
//    call gravity_bicgstab_start_2
//
// --------------------
// gravity_bicgstab_start_2()
// --------------------
//
//    if (is_singular)
//      B = B - bs/bc;
//
//    R = B
//    R0 = R
//    P = R
//    beta_d_ = DOT(R, R) ==> r_gravity_bicgstab_start_3
//
// --------------------
// r_gravity_bicgstab_start_3()
// --------------------
//
//    receive beta_d_
// 
//    call gravity_bicgstab_loop_0()
//
// ==================================================
//
// --------------------
// gravity_bicgstab_loop_0()
// --------------------
//
//    if (iter_ == 0) {
//       rho0_ = sqrt(beta_d_)
//       if (rho0_ == 0.0)  rho0_ = 1.0
//       err_ = sqrt(beta_d_) / rho0_
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
//      ==> gravity_bicgstab_end(return_converged)
//
//    if (iter_ >= iter_max_) {
//       ==> gravity_bicgstab_end(return_error_max_iter_reached);
//    }
// 
//    refresh(P) ==> r_gravity_bicgstab_loop_1
//
// --------------------
// r_gravity_bicgstab_loop_1()
// --------------------
//
//    call gravity_bicgstab_loop_2
//
// --------------------
// gravity_bicgstab_loop_2()
// --------------------
//
//    Y = SOLVE(M,P)
//
//    refresh(Y) ==> r_gravity_bicgstab_loop_2
//
// --------------------
// r_gravity_bicgstab_loop_3()
// --------------------
//  
//    call gravity_bicgstab_loop_4
// 
// --------------------
// gravity_bicgstab_loop_4()
// --------------------
//
//    V = MATVEC(A,Y)
//
//    vr0_ = DOT(V, R0) ==> r_gravity_bicgstab_loop_5
//
// --------------------
// r_gravity_bicgstab_loop_5()
// --------------------
//
//    receive vr0_
//
//    call gravity_bicgstab_loop_6
//
// --------------------
// gravity_bicgstab_loop_6()
// --------------------
//
//    alpha_ = beta_d_ / ( vr0_ );
//    Q = R - alpha_*V
//    X = X + alpha_*Y
//
//    refresh(Q) ==> r_gravity_bicgstab_loop_7
//
// --------------------
// r_gravity_bicgstab_loop_7()
// --------------------
//  
//    call gravity_bicgstab_loop_8
//
// --------------------
// gravity_bicgstab_loop_8()
// --------------------
//  
//    Y = SOLVE(M,Q)
//
//    refresh(Y) ==> r_gravity_bicgstab_loop_9
// 
// --------------------
// r_gravity_bicgstab_loop_9()
// --------------------
//  
//    call gravity_bicgstab_loop_10
// 
// --------------------
// gravity_bicgstab_loop_10()
// --------------------
//
//    V = MATVEC(A,Y)
//
//    omega_d_ = DOT(U, U) ==> r_gravity_bicgstab_loop_11
//    omega_n_ = DOT(U, Q) ==> r_gravity_bicgstab_loop_11
//
// --------------------
// r_gravity_bicgstab_loop_11()
// --------------------
//
//    receive omega_d_ and omega_n_
//
//    call gravity_bicgstab_loop_12
//
// --------------------
// gravity_bicgstab_loop_12()
// --------------------
//
//    if (omega_d_ == 0.0)  omega_d_ = 1.0
//    if ( omega_n_ == 0.0 ) 
//       gravity_bicgstab_end(return_error_omega_eq_0)
//
//    omega_ = omega_n_ / omega_d_
//    X = X + omega_*Y
//    R = Q - omega_*U
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
//    call gravity_bicgstab_loop_14
//
// --------------------
// gravity_bicgstab_loop_14()
// --------------------
//
//    if ( beta_n_ == 0.0 )
//       gravity_bicgstab_end(return_error_beta_n_eq_0)
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
//    call gravity_bicgstab_loop_0
//
//
// ==================================================
//
// --------------------
// gravity_bicgstab_end(return_)
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
(const FieldDescr* field_descr, int rank,
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
    bs_(0.0), bc_(0.0),
    id_refresh_matvec_(-1)
{

  M_ = (diag_precon) ? (Matrix*)(new EnzoMatrixDiagonal) 
                     : (Matrix*)(new EnzoMatrixIdentity);

  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  ix_ = field_descr->field_id("X");
  ir_ = field_descr->field_id("R");
  ir0_ = field_descr->field_id("R0");
  ip_ = field_descr->field_id("P");
  iy_ = field_descr->field_id("Y");
  iv_ = field_descr->field_id("V");
  iq_ = field_descr->field_id("Q");
  iu_ = field_descr->field_id("U");

  // Initialize default Refresh
  const int num_fields = field_descr->field_count();

  const int ir = add_refresh(1, rank-1, neighbor_leaf, sync_barrier);
  //  refresh(ir)->add_field(idensity_);
  refresh(ir)->add_all_fields(num_fields);

  // Initialize matvec Refresh
  id_refresh_matvec_ = add_refresh(1, rank-1, neighbor_leaf, sync_barrier);
  refresh(id_refresh_matvec_)->add_all_fields(num_fields);
  //  refresh(id_refresh_matvec_)->add_field(ir_);

}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::compute(Block* block) throw() {
  
  Field field = block->data()->field();

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (block);

  field.size(&nx_, &ny_, &nz_);
  field.dimensions(idensity_, &mx_, &my_, &mz_);
  field.ghost_depth(idensity_, &gx_, &gy_, &gz_);

  int precision = field.precision(idensity_);

  if      (precision == precision_single)    compute_<float>      (enzo_block);
  else if (precision == precision_double)    compute_<double>     (enzo_block);
  else if (precision == precision_quadruple) compute_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
}

//======================================================================

extern CkReduction::reducerType r_method_gravity_bicgstab_type;

template<class T> void EnzoMethodGravityBiCGStab::compute_(EnzoBlock* enzo_block) throw() {

  // initialize iter to 0
  iter_ = 0;

  // construct B <right-hand side>
  // initialize X <initial solution X0> to zero
  Data* data  = enzo_block->data();
  Field field = data->field();

  if (enzo_block->is_leaf()) {

    T* density = (T*) field.values(idensity_);
    T* B       = (T*) field.values(ib_);
    T* X       = (T*) field.values(ix_);
    T* R       = (T*) field.values(ir_);

    //   X = 0   [Q: is this necessary?  couldn't we reuse the solution from the previous solve?]
    //   B = -h^2 * 4 * PI * G * density
    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  X[i] = 0.0;
	  B[i] = - 4.0 * (cello::pi) * grav_const_ * density[i];
	}
      }
    }
  }

  // set bs_ = SUM(B)   ==> r_gravity_bicgstab_start_1
  // set bc_ = COUNT(B) ==> r_gravity_bicgstab_start_1
  long double reduce[2] = {0.0, 0.0};

  if (enzo_block->is_leaf()) {
    reduce[0] = sum_(B);
    reduce[1] = count_(B);
  }

  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_start_1<T>(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute(2*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_start_1(CkReductionMsg* msg) {

  // EnzoBlock accumulate global contributions to SUM(B) and COUNT(B)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());
  long double* data = (long double*) msg->getData();
  method->set_bs( data[0] );
  method->set_bc( data[1] );
  delete msg;

  // call gravity_bicgstab_start_2 to continue
  method->gravity_bicgstab_start_2(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_start_2(EnzoBlock* enzo_block) throw() {

  Data* data = enzo_block->data();
  Field field = data->field();

  // if (is_singular)   B = B - bs/bc;
  //
  // R = B
  // R0 = R
  // P = R
  if (enzo_block->is_leaf()) {

    T* B       = (T*) field.values(ib_);
    T* R0      = (T*) field.values(ir0_);
    T* P       = (T*) field.values(ip_);
    T* R       = (T*) field.values(ir_);

    // shift B if the problem is singular
    if (is_singular_) {

      T shift = -bs_ / bc_;

      for (int iz=0; iz<mz_; iz++) {
	for (int iy=0; iy<my_; iy++) {
	  for (int ix=0; ix<mx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    B[i] = shift + B[i];
	  }
	}
      }
    }

    // set R = R0 = P = B
    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  R[i] = R0[i] = P[i] = B[i];
	}
      }
    }

  }

  // beta_d_ = DOT(R, R) ==> r_gravity_bicgstab_start_3
  long double reduce = 0.0;
  if (enzo_block->is_leaf())
    reduce = dot_(R,R);
  
  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_start_3<T>(NULL), 
		      this->proxy_array());

  this->contribute(sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_start_3(CkReductionMsg* msg) {

  // EnzoBlock accumulate global contributions to DOT(B,B)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  long double* data = (long double*) msg->getData();
  method->set_beta_d( data[0] );
  delete msg;

  // call gravity_bicgstab_loop_0 to begin solver loop
  method->gravity_bicgstab_loop_0<T>(this);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_0(EnzoBlock* enzo_block) throw() {

  cello::check(beta_d_,"beta_d_",__FILE__,__LINE__);

  // initialize/update current error, store error statistics
  if (iter_ == 0) {
    rho0_ = sqrt(beta_d_)
    if (rho0_ == 0.0)  rho0_ = 1.0
    err_ = sqrt(beta_d_) / rho0_;
    err_min_ = err_;
    err_max_ = err_;
  } else {
    err_ = sqrt(rr_) / rho0_;
    err_min_ = std::min(err_, err_min_);
    err_max_ = std::max(err_, err_max_);
  }

  // output solution progress (iteration, residual, etc)
  if (enzo_block->index().is_root()) {

    Monitor* monitor = enzo_block->simulation()->monitor();

    if (iter_ == 0) {
      monitor->print("Enzo", "BiCGStab iter %04d  rho0 %g",
		     iter_,(double)(rho0_));
    }

    if (monitor_iter_ && (iter_ % monitor_iter_) == 0 ) {
      monitor_output_ (enzo_block);
    }
  }

  // check for convergence
  if (err_ < res_tol_) {

    gravity_bicgstab_end<T>(enzo_block, return_converged);


  // check for failure
  } else if (iter_ >= iter_max_)  {

    gravity_bicgstab_end<T>(enzo_block, return_error_max_iter_reached);


  // refresh P and call r_gravity_bicgstab_loop_1 to handle re-entry
  } else {

    this->refresh(1)->set_active(enzo_block->is_leaf());
    enzo_block->refresh_enter(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_1(NULL),
			      this->refresh(1));

  }
}

//----------------------------------------------------------------------

void EnzoBlock::r_gravity_bicgstab_loop_1() {
  /// re-entry into gravity_bicgstab_loop_2

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  if      (precision == precision_single)    
    method->gravity_bicgstab_loop_2<float>(enzo_block);
  else if (precision == precision_double)    
    method->gravity_bicgstab_loop_2<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->gravity_bicgstab_loop_2<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_2(EnzoBlock* enzo_block) throw() {

  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;
 
  // Y = SOLVE(M, P)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);
    M_->matvec(iy_, ip_, enzo_block);
  }

  // refresh Y and call r_gravity_bicgstab_loop_3 to handle re-entry
  this->refresh(1)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_enter(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_3(NULL),
			    this->refresh(1));

}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_3(CkReductionMsg* msg) {
  /// re-entry into gravity_bicgstab_loop_4

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  if      (precision == precision_single)    
    method->gravity_bicgstab_loop_4<float>(enzo_block);
  else if (precision == precision_double)    
    method->gravity_bicgstab_loop_4<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->gravity_bicgstab_loop_4<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_4(EnzoBlock* enzo_block) throw() {

  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;
 
  // V = MATVEC(A,Y)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);
    A_->matvec(iv_, iy_, enzo_block);
  }

  // vr0_ = DOT(V, R0) ==> r_gravity_bicgstab_loop_5
  long double reduce = 0.0;
  if (enzo_block->is_leaf()) {
    T* R0 = (T*) field.values(ir0_);
    T* V  = (T*) field.values(iv_);
    reduce = dot_(V,R0);
  }
  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_5<T>(NULL), 
		enzo_block->proxy_array());

  enzo_block->contribute(sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
  
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_5(CkReductionMsg* msg) {

  // EnzoBlock accumulate global contributions to DOT(V,R0)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  method->set_vr0( ((long double*)msg->getData())[0] );
  delete msg;

  // call gravity_bicgstab_loop_6 to continue
  method->gravity_bicgstab_loop_6<T>(this);
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_6(EnzoBlock* enzo_block) throw() {

  cello::check(vr0_,"vr0_",__FILE__,__LINE__);

  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;

  if (enzo_block->is_leaf()) {

    data->field_cell_width(&hx,&hy,&hz);
    T* Q = (T*) field.values(iq_);
    T* R = (T*) field.values(ir_);
    T* V = (T*) field.values(iv_);
    T* X = (T*) field.values(ix_);
    T* Y = (T*) field.values(iy_);

    // alpha_ = beta_d_ / vr0_;
    alpha_ = beta_d_ / vr0_;

    // Q = -alpha_*V + R
    zaxpy_(Q, -alpha_, V, R);

    // X = alpha_*Y + X
    zaxpy_(X, alpha_, Y, X);

  }

  // refresh Q and call r_gravity_bicgstab_loop_7 to handle re-entry
  this->refresh(1)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_enter(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_7(NULL),
			    this->refresh(1));

}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_7(CkReductionMsg* msg) {
  /// re-entry into gravity_bicgstab_loop_8

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  if      (precision == precision_single)    
    method->gravity_bicgstab_loop_8<float>(enzo_block);
  else if (precision == precision_double)    
    method->gravity_bicgstab_loop_8<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->gravity_bicgstab_loop_8<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_8(EnzoBlock* enzo_block) throw() {

  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;
 
  // Y = SOLVE(M, Q)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);
    M_->matvec(iy_, iq_, enzo_block);
  }

  // refresh Y and call r_gravity_bicgstab_loop_9 to handle re-entry
  this->refresh(1)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_enter(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_9(NULL),
			    this->refresh(1));

}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_9(CkReductionMsg* msg) {
  /// re-entry into gravity_bicgstab_loop_10

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  if      (precision == precision_single)    
    method->gravity_bicgstab_loop_10<float>(enzo_block);
  else if (precision == precision_double)    
    method->gravity_bicgstab_loop_10<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->gravity_bicgstab_loop_10<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_10(EnzoBlock* enzo_block) throw() {

  Data* data = enzo_block->data();
  Field field = data->field();
  double hx,hy,hz;

  // V = MATVEC(A,Y)
  if (enzo_block->is_leaf()) {
    data->field_cell_width(&hx,&hy,&hz);
    A_->matvec(iv_, iy_, enzo_block);
  }

  // omega_d_ = DOT(U, U) ==> r_gravity_bicgstab_loop_11
  // omega_n_ = DOT(U, Q) ==> r_gravity_bicgstab_loop_11
  long double reduce[2] = {0.0, 0.0};
  if (enzo_block->is_leaf()) {
    T* U  = (T*) field.values(iu_);
    T* Q  = (T*) field.values(iq_);
    reduce[0] = dot_(U,U);
    reduce[1] = dot_(U,Q);
  }

  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_11<T>(NULL), 
		enzo_block->proxy_array());

  enzo_block->contribute(2*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
    
}


//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_11(CkReductionMsg* msg) {

  // EnzoBlock accumulate global contributions to DOT(U,U) and DOT(U,Q)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  long double* data = (long double*) msg->getData();
  method->set_omega_d( data[0] );
  method->set_omega_n( data[1] );

  delete msg;

  // call gravity_bicgstab_loop_12 to continue
  method->gravity_bicgstab_loop_12<T>(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_12(EnzoBlock* enzo_block) throw() {

  cello::check(omega_d_,"omega_d_",__FILE__,__LINE__);
  cello::check(omega_n_,"omega_n_",__FILE__,__LINE__);

  Data* data = enzo_block->data();
  Field field = data->field();

  // fix omega_d_ if needed
  if (omega_d_ == 0.0)  omega_d_ = 1.0;
  
  // compute omega_ = omega_n_ / omega_d_
  omega_ = omega_n_ / omega_d_;

  // check for error
  if ( omega_ == 0.0 )
    gravity_bicgstab_end<T>(enzo_block, return_error_omega_eq_0);


  // update: X = X + omega_*Y
  //         R = Q - omega_*U
  if (enzo_block->is_leaf()) {

    T* X = (T*) field.values(ix_);
    T* Y = (T*) field.values(iy_);
    T* R = (T*) field.values(ir_);
    T* Q = (T*) field.values(iq_);
    T* U = (T*) field.values(iu_);
    
    // X = omega_*Y + X
    zaxpy_(X, omega_, Y, X);

    // R = -omega_*U + Q
    zaxpy_(R, -omega_, U, Q);

  }

  // rr_     = DOT(R, R)  ==> r_gravity_bicgstab_loop_13
  // beta_n_ = DOT(R, R0) ==> r_gravity_bicgstab_loop_13
  long double reduce[2] = {0.0, 0.0};
  if (enzo_block->is_leaf()) {
    T* R  = (T*) field.values(ir_);
    T* R0 = (T*) field.values(ir0_);
    reduce[0] = dot_(R,R);
    reduce[1] = dot_(R,R0);
  }

  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_13<T>(NULL), 
		enzo_block->proxy_array());
  
  enzo_block->contribute(2*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
    
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_13(CkReductionMsg* msg) {

  // EnzoBlock accumulate global contributions to DOT(R,R) and DOT(R,R0)
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  long double* data = (long double*) msg->getData();

  method->set_rr(     data[0] );
  method->set_beta_n( data[1] );

  delete msg;

  // call gravity_bicgstab_loop_14 to continue
  method->gravity_bicgstab_loop_14<T>(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_loop_14(EnzoBlock* enzo_block) throw() {

  cello::check(rr_,    "rr_",    __FILE__,__LINE__);
  cello::check(beta_n_,"beta_n_",__FILE__,__LINE__);

  Data* data = enzo_block->data();
  Field field = data->field();

  // check for failure via beta_n_ == 0
  if (beta_n_ == 0.0)
    gravity_bicgstab_end<T>(enzo_block, return_error_beta_n_eq_0);
  
  // update: beta_ = (beta_n_/beta_d_)*(alpha_/omega_)
  beta_ = (beta_n_/beta_d_)*(alpha_/omega_);

  // update: P = R + beta_*( P - omega_*V )
  if (enzo_block->is_leaf()) {

    T* P = (T*) field.values(ip_);
    T* R = (T*) field.values(ir_);
    T* V = (T*) field.values(iv_);

    // P = beta_*P + R
    zaxpy_(P, beta_, P, R);

    // P = beta_*omega_*V + P
    zaxpy_(P, beta_*omega_, V, P);

  }

  // contribute to iteration counter ==> r_gravity_bicgstab_loop_15
  int iter = iter_ + 1;

  CkCallback callback(CkIndex_EnzoBlock::r_gravity_bicgstab_loop_15<T>(NULL), 
		      enzo_block->proxy_array());
    
  enzo_block->contribute(sizeof(int), &iter, 
			 CkReduction::max_int, callback);

}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_gravity_bicgstab_loop_15(CkReductionMsg* msg) {

  // EnzoBlock accumulate global contributions to iter
  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  method->set_iter ( ((int*)msg->getData())[0] );

  delete msg;

  // call gravity_bicgstab_loop_0 to continue to next iteration
  method->gravity_bicgstab_loop_0<T>(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::gravity_bicgstab_end(EnzoBlock* enzo_block, int retval) throw () {
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
  Data* data = enzo_block->data();
  Field field = data->field();

  if (enzo_block->is_leaf()) {

    T* X         = (T*) field.values(ix_);
    T* potential = (T*) field.values(ipotential_);

    if (enzo_block->index().is_root()) 
      monitor_output_(enzo_block);

    copy_(potential, X, mx_, my_, mz_, is_leaf_);

    bool symmetric;
    int order;
    EnzoComputeAcceleration compute_acceleration(field.field_descr(),
						 rank_, symmetric=true,
						 order=2);
    compute_acceleration.compute(enzo_block);
  }

  enzo_block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::monitor_output_(EnzoBlock* enzo_block) throw() {

  Monitor* monitor = enzo_block->simulation()->monitor();

  monitor->print("Enzo", "BiCGStab iter %04d  rr %g [%g %g]",
		 iter_,
		 (double)(rr_    / rho0_),
		 (double)(rr_min_/ rho0_),
		 (double)(rr_max_/ rho0_));
}

//======================================================================

template<class T> long double EnzoMethodGravityBiCGStab::dot_(const T* X, const T* Y) const throw() {

  const int i0 = gx_ + mx_*(gy_ + my_*gz_);

  long double value = 0.0;
  for (int iz=0; iz<nz_; iz++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + (ix + mx_*(iy + my_*iz));
	value += X[i]*Y[i];
      }
    }
  }

  return value;
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::zaxpy_(T* Z, double a, const T* X, const T* Y) const throw() {

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	Z[i] = a * X[i] + Y[i];
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T> long double EnzoMethodGravityBiCGStab::sum_(const T* X) const throw() {

  const int i0 = gx_ + mx_*(gy_ + my_*gz_);

  long double value = 0.0;
  for (int iz=0; iz<nz_; iz++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + (ix + mx_*(iy + my_*iz));
	value += X[i];
      }
    }
  }

  return value;
}

//----------------------------------------------------------------------

long double EnzoMethodGravityBiCGStab::count_() const throw() {
  return 1.0*nx_*ny_*nz_;
}

//----------------------------------------------------------------------
