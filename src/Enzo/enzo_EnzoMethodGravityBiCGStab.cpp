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
// function [x, error, iter, flag] = bcgs_opt(A, x, b, M, max_it, tol)
// % Usage: [x, error, iter, flag] = bcgs_opt(A, x, b, M, max_it, tol)
// %
// % Solves the linear system Ax=b using the BiConjugate Gradient
// % Stabilized Method with right preconditioning.  We assume an
// % initial guess of x=0.
// %
// % input   A        REAL matrix
// %         b        REAL right hand side vector
// %         M        REAL preconditioner matrix
// %         max_it   INTEGER maximum number of iterations
// %         tol      REAL error tolerance
// %
// % output  x        REAL solution vector
// %         error    REAL error norm
// %         iter     INTEGER number of iterations performed
// %         flag     INTEGER: 0 = solution found to tolerance
// %                           1 = no convergence given max_it
// %                          -1 = breakdown: rho = 0
// %                          -2 = breakdown: omega = 0
//
//   % initialize outputs
//   iter = 0;
//   flag = 0;
//
//   % compute initial residual, initialize temporary vectors
//   r = b - A*x;                       % pt-to-pt communication
//   rstar = r;
//   p = r;
//   
//   % inner-products at initialization: <b,b> and <r,r>
//   rho0 = dot(b,b);                   % global reduction
//   beta_d = dot(r,r);                 % global reduction
//   
//   % store ||b|| for relative error checks
//   rho0 = sqrt(rho0);
//   if (rho0 == 0.0), rho0 = 1.0; end
//
//   % check initial relative residual to see if work is finished
//   error = sqrt(beta_d)/rho0;
//   if (error < tol), return; end
//
//   % begin iteration
//   for iter = 1:max_it
//
//      % first half of step
//      y = M\p;                        % ? communication
//      v = A*y;                        % pt-to-pt communication
//      alpha = beta_d / dot(v,rstar);  % global reduction
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
//      beta_n = dot(r,rstar);          % global reduction
//      
//      % check for failure/convergence
//      if (beta_n == 0), flag = -1; return; end
//      error = sqrt(rr)/rho0;
//      if (error < tol), return; end
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
// BiCGStab partitioned along parallel communication / synchronization steps
//
// --------------------
// bicgstab_begin()
// --------------------
//
//    return_ = return_unknown;
// 
//    B = <right-hand side>
//    X = <initial solution X0>
//    R = MATVEC(A,x)  ==> bicgstab_start_1
//
// --------------------
// bicgstab_start_1()
// --------------------
//
//    Rstar = R
//    P = R
//    rho0_   = DOT(B, B) ==> bicgstab_start_2
//    beta_d_ = DOT(R, R) ==> bicgstab_start_2
//
// --------------------
// bicgstab_start_2()
// --------------------
//
//    rho0 = sqrt(rho0_)
//    if (rho0 == 0.0)  rho0 = 1.0
//    error = sqrt(beta_d_) / rho0;
//
//    if (error < tol) {
//      ==> bicgstab_loop_end(return_converged_)
//    }
//    
//    iter = 1
//    
//    ==> bicgstab_loop_begin()
//
// ==================================================
//
// --------------------
// bicgstab_loop_begin()
// --------------------
//
//    if (iter >= max_it) {
//       ==> bicgstab_loop_end(return_error_not_converged_);
//    }
// 
//    Y = SOLVE (M, P) ==> bicgstab_loop_1
// 
// --------------------
// bicgstab_loop_1()
// --------------------
//  
//    V = MATVEC (A,Y) ==> bicgstab_loop_2
// 
// --------------------
// bicgstab_loop_2()
// --------------------
//
//    vrs_ = DOT(V, Rstar) ==> bicgstab_loop_3
//
// --------------------
// bicgstab_loop_3()
// --------------------
//
//    alpha_ = beta_d_ / ( vrs_ );
//    Q = R - alpha_*V
//    X = X + alpha_*Y
//
//    Y = SOLVE (M, Q) ==> bicgstab_loop_4
//
// --------------------
// bicgstab_loop_4()
// --------------------
//
//    U = MATVEC (A,Y) ==> bicgstab_loop_5
// 
// --------------------
// bicgstab_loop_5()
// --------------------
//
//    omega_d_ = DOT(U, U) ==> bicgstab_loop_6
//    omega_n_ = DOT(U, Q) ==> bicgstab_loop_6
//
// --------------------
// bicgstab_loop_6()
// --------------------
//
//    if (omega_d_ == 0.0)  omega_d_ = 1.0
//    omega_ = omega_n_ / omega_d_
//    if ( omega_ == 0.0 ) {
//       bicgstab_loop_end(return_error_omega_eq_0)
//    }
// 
//    X = X + omega_*Y
//    R = Q - omega_*U
//
//    rr_     = DOT(R, R)     ==> bicgstab_loop_7
//    beta_n_ = DOT(R, Rstar) ==> bicgstab_loop_7
//
// --------------------
// bicgstab_loop_7()
// --------------------
//
//    if ( beta_n_ == 0.0 ) {
//       bicgstab_loop_end(return_error_beta_n_eq_0)
//    }
//    error = sqrt(rr_) / rho0;
//
//    if ( error <= tol ) {
//       bicgstab_loop_end(return_converged_)
//    }
//
//    beta = (beta_n_/beta_d_)*(alpha_/omega_)
//    P = R + beta*( P - omega_*v )
//    iter = iter + 1
//
//    ==> bicgstab_loop_begin()
//
// ==================================================
//
// --------------------
// bicgstab_loop_end(return_)
// --------------------
//
// if (return_ == return_converged) {
//    NEXT()
// } else {
//    ERROR (return_)
// }
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
    rr0_(0),
    rr_min_(0),rr_max_(0),
    idensity_(0),  ipotential_(0),
    ib_(0), ix_(0), ir_(0), id_(0), iy_(0), iz_(0),
    is_leaf_(false),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    rr_(0.0), rz_(0.0), rz2_(0.0), dy_(0.0), bs_(0.0), bc_(0.0),
    id_refresh_matvec_(-1)
{

  M_ = (diag_precon) ? (Matrix*)(new EnzoMatrixDiagonal) 
    :                  (Matrix*)(new EnzoMatrixIdentity);

  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  id_ = field_descr->field_id("D");
  ir_ = field_descr->field_id("R");
  ix_ = field_descr->field_id("X");
  iy_ = field_descr->field_id("Y");
  iz_ = field_descr->field_id("Z");

  // Initialize default Refresh
  const int ir = add_refresh(1,rank-1,sync_barrier);
  refresh(ir)->add_field(idensity_);

  // Initialize matvec Refresh
  id_refresh_matvec_ = add_refresh(1,rank-1,sync_barrier);
  refresh(id_refresh_matvec_)->add_field(ir_);

}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::compute(Block* block) throw() {
  
  set_leaf(block);

  Field field = block->data()->field();

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (block);

  field.size                (&nx_,&ny_,&nz_);
  field.dimensions(idensity_,&mx_,&my_,&mz_);
  field.ghost_depth    (idensity_,&gx_,&gy_,&gz_);

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
  //     X = initial guess
  //     B = right-hand side
  //     R = B - A*X
  //     solve(M*Z = R)
  //     D = Z
  //     shift (B)

  set_leaf(enzo_block);

  iter_ = 0;

  Data* data = enzo_block->data();
  Field field = data->field();

  T* density = (T*) field.values(idensity_);
    
  T* B = (T*) field.values(ib_);
  T* X = (T*) field.values(ix_);
  T* R = (T*) field.values(ir_);

  if (is_leaf_) {

    //   - X = 0
    //   - B = -h^2 * 4 * PI * G * density
    //   - R = P = B ( residual with X = 0);
    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  X[i] = 0.0;
	  B[i] = - 4.0 * (cello::pi) * grav_const_ * density[i];
	  R[i] = B[i];
	}
      }
    }
  }

  M_->matvec(id_,ir_,enzo_block);
  M_->matvec(iz_,ir_,enzo_block);

  long double reduce[3];

  reduce[0] = dot_(R,R);
  reduce[1] = sum_(B);
  reduce[2] = count_(B);

  CkCallback callback(CkIndex_EnzoBlock::r_bicgstab_loop_0a<T>(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute(3*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_bicgstab_loop_0a(CkReductionMsg* msg) {
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> refresh P for AP = MATVEC (A,P)

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  long double* data = (long double*) msg->getData();

  method->set_rr( data[0] );
  method->set_bs( data[1] );
  method->set_bc( data[2] );

  delete msg;
  
  // Refresh if Block is a leaf
  method->refresh(1)->set_active(is_leaf());
  refresh_enter(CkIndex_EnzoBlock::r_enzo_matvec(NULL),
		method->refresh(1));
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_bicgstab_loop_0b(CkReductionMsg* msg) {
/// ==> refresh P for AP = MATVEC (A,P)

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  method->set_iter( ((int*)msg->getData())[0] );

  // Refresh if Block is a leaf
  method->refresh(1)->set_active(is_leaf());
  refresh_enter(CkIndex_EnzoBlock::r_enzo_matvec(NULL),
		method->refresh(1));
}

//----------------------------------------------------------------------

// void EnzoBlock::enzo_matvec_() {

//   EnzoMethodGravityBiCGStab* method = 
//     static_cast<EnzoMethodGravityBiCGStab*> (this->method());

//   Field field = data()->field();
//   int precision = field.precision(field.field_id("density")); // assuming 

//   EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

//   if (precision == precision_single)    
//     method->bicgstab_shift_1<float>(enzo_block);
//   else if (precision == precision_double)    
//     method->bicgstab_shift_1<double>(enzo_block);
//   else if (precision == precision_quadruple) 
//     method->bicgstab_shift_1<long double>(enzo_block);
//   else 
//     ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", precision);
// }

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::bicgstab_shift_1(EnzoBlock* enzo_block) throw() {

  cello::check(rr_,"rr_",__FILE__,__LINE__);
  cello::check(bs_,"bs_",__FILE__,__LINE__);
  cello::check(bc_,"bc_",__FILE__,__LINE__);

  set_leaf(enzo_block);

  Data* data = enzo_block->data();
  Field field = data->field();

  T* B  = (T*) field.values(ib_);
  T* R  = (T*) field.values(ir_);

  if (iter_ == 0 && is_singular_)  {

    // shift rhs B by projection of B onto e: B~ <== B - (e*eT)/(eT*e) b
    // eT*e == n === zone count (bc)
    // eT*b == sum_i=1,n B[i]

    T shift = -bs_ / bc_;
    shift_ (R,shift,R);
    shift_ (B,shift,B);

    M_->matvec(id_,ir_,enzo_block);
    M_->matvec(iz_,ir_,enzo_block);
  } 

  long double reduce;

  reduce = dot_(R,R);

  CkCallback callback(CkIndex_EnzoBlock::r_bicgstab_shift_1<T>(NULL), 
		enzo_block->proxy_array());

  enzo_block->contribute(sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
    
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_bicgstab_shift_1(CkReductionMsg* msg) {

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  method->set_rr( ((long double*)msg->getData())[0] );

  delete msg;

  method->bicgstab_loop_2<T>(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::bicgstab_loop_2(EnzoBlock* enzo_block) throw() {

  set_leaf(enzo_block);

  cello::check(rr_,"rr_",__FILE__,__LINE__);

  if (iter_ == 0) {
    rr0_ = rr_;
    rr_min_ = rr_;
    rr_max_ = rr_;
  } else {
    rr_min_ = std::min(rr_min_,rr_);
    rr_max_ = std::max(rr_max_,rr_);
  }

  if (enzo_block->index().is_root()) {

    Monitor* monitor = enzo_block->simulation()->monitor();

    if (iter_ == 0) {
      monitor->print("Enzo", "BiCGStab iter %04d  rr0 %g",
		     iter_,(double)(rr0_));
    }

    if (monitor_iter_ && (iter_ % monitor_iter_) == 0 ) {
      monitor_output_ (enzo_block);
    }
  }

  if (rr_ / rr0_ < res_tol_) {

    bicgstab_end<T>(enzo_block,return_converged);

  } else if (iter_ >= iter_max_)  {

    bicgstab_end<T>(enzo_block,return_error_max_iter_reached);

  } else {

    Data* data = enzo_block->data();
    Field field = data->field();

    T* D = (T*) field.values(id_);
    T* Y = (T*) field.values(iy_);
    T* R = (T*) field.values(ir_);
    T* Z = (T*) field.values(iz_);

    double hx,hy,hz;
    data->field_cell_width(&hx,&hy,&hz);

    A_->matvec(iy_,id_,enzo_block);

    long double reduce[3];

    reduce[0] = dot_(R,R);
    reduce[1] = dot_(R,Z);
    reduce[2] = dot_(D,Y);

    CkCallback callback(CkIndex_EnzoBlock::r_bicgstab_loop_3<T>(NULL), 
		  enzo_block->proxy_array());

    enzo_block->contribute(3*sizeof(long double), &reduce, 
			   r_method_gravity_bicgstab_type, callback);
  }
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_bicgstab_loop_3(CkReductionMsg* msg) {

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  long double* data = (long double*) msg->getData();

  method->set_rr(data[0]);
  method->set_rz(data[1]);
  method->set_dy(data[2]);
  
  delete msg;

  method->bicgstab_loop_4<T>(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::bicgstab_loop_4(EnzoBlock* enzo_block) throw () {
  //  a = rz / dy;
  //  X = X + a*D;
  //  R = R - a*Y;
  //  solve(M*Z = R)
  //  rz2 = dot(R,Z)
  //  b = rz2 / rz;
  //  D = Z + b*D;
  //  rz = rz2;

  set_leaf(enzo_block);

  cello::check(rr_,"rr_",__FILE__,__LINE__);
  cello::check(rz_,"rz_",__FILE__,__LINE__);
  cello::check(dy_,"dy_",__FILE__,__LINE__);

  Data* data = enzo_block->data();
  Field field = data->field();

  T* X = (T*) field.values(ix_);
  T* D = (T*) field.values(id_);
  T* R = (T*) field.values(ir_);
  T* Y = (T*) field.values(iy_);
  T* Z = (T*) field.values(iz_);

  T a = rz_ / dy_;

  cello::check(a,"a",__FILE__,__LINE__);

  zaxpy_(X,  a ,D,X);
  zaxpy_(R, -a, Y,R);

  M_->matvec(iz_,ir_,enzo_block);

  long double reduce[3];

  reduce[0] = dot_(R,Z);
  reduce[1] = sum_(R);
  reduce[2] = sum_(X);

  CkCallback callback(CkIndex_EnzoBlock::r_bicgstab_loop_5<T>(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute(3*sizeof(long double), &reduce, 
			 r_method_gravity_bicgstab_type, callback);
}

//----------------------------------------------------------------------

template<class T> void EnzoBlock::r_bicgstab_loop_5(CkReductionMsg* msg) {
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> bicgstab_loop_6

  EnzoMethodGravityBiCGStab* method = 
    static_cast<EnzoMethodGravityBiCGStab*> (this->method());

  long double* data = (long double*) msg->getData();

  method->set_rz2(data[0]);
  method->set_rs (data[1]);
  method->set_xs (data[2]);

  delete msg;

  method->bicgstab_loop_6<T>(this);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::bicgstab_loop_6(EnzoBlock* enzo_block) throw () {
  //  rz2 = dot(R,Z)
  //  b = rz2 / rz;
  //  D = Z + b*D;
  //  rz = rz2;

  set_leaf(enzo_block);

  cello::check(rz2_,"rz2_",__FILE__,__LINE__);
  cello::check(rs_,"rs_",__FILE__,__LINE__);
  cello::check(xs_,"xs_",__FILE__,__LINE__);

  Field field = enzo_block->data()->field();

  if (is_singular_)  {

    // shift rhs B by projection of B onto e: B~ <== B - (e*eT)/(eT*e) b
    // eT*e == n === zone count (bc)
    // eT*b == sum_i=1,n B[i]

    T* X = (T*) field.values(ix_);
    T* R = (T*) field.values(ir_);
    shift_(X, T(-xs_/bc_), X);
    shift_(R, T(-rs_/bc_), R);
  }

  T* D = (T*) field.values(id_);
  T* Z = (T*) field.values(iz_);

  T b = rz2_ / rz_;

  cello::check(b,"b",__FILE__,__LINE__);

  zaxpy_(D, b, D, Z);

  int iter = iter_ + 1;

  CkCallback callback(CkIndex_EnzoBlock::r_bicgstab_loop_0b<T>(NULL), 
		      enzo_block->proxy_array());
    
  enzo_block->contribute(sizeof(int), &iter, 
			 CkReduction::max_int, callback);

}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::bicgstab_end(EnzoBlock* enzo_block, int retval) throw () {
  ///    if (return == return_converged) {
  ///       potential = X
  ///       ==> bicgstab_exit()
  ///    } else {
  ///       ERROR (return-)
  ///    }

  set_leaf(enzo_block);

  Data* data = enzo_block->data();
  Field field = data->field();

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
  enzo_block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::monitor_output_(EnzoBlock* enzo_block) throw() {

  Monitor* monitor = enzo_block->simulation()->monitor();

  monitor->print("Enzo", "BiCGStab iter %04d  rr %g [%g %g]",
		 iter_,
		 (double)(rr_    / rr0_),
		 (double)(rr_min_/ rr0_),
		 (double)(rr_max_/ rr0_));
}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::bicgstab_exit_() throw() {
  // deallocate temporary vectors
}

//======================================================================

template<class T> long double EnzoMethodGravityBiCGStab::dot_(const T* X, const T* Y) const throw() {

  if (!is_leaf_) return 0.0;

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

  if (!is_leaf_) return;

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

  if (!is_leaf_) return 0.0;

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

template<class T> void EnzoMethodGravityBiCGStab::shift_(T* X, const T a, const T* Y) const throw() {

  if (!is_leaf_) return;

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] = a + Y[i];
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T> void EnzoMethodGravityBiCGStab::scale_(T* Y, T a, const T* X) const throw() {

  if (!is_leaf_) return;

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	Y[i] = a * X[i];
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T> int EnzoMethodGravityBiCGStab::count_(T* X) const throw() {
  return (is_leaf_) ? nx_*ny_*nz_ : 0;
}

//----------------------------------------------------------------------
