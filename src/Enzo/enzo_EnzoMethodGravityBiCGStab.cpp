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
//   beta_n = beta_d;
//   error = sqrt(beta_d)/rho0;
//   if (error < tol), return; end
//
//   % begin iteration
//   for iter = 1:max_it
//
//      % first half of Bi-CG step
//      y = M\p;
//      v = A*y;                        % pt-to-pt communication
//      alpha = beta_n / dot(v,rstar);  % global reduction
//      q = r - alpha*v;
//      x = x + alpha*y;
//
//      % stabilization portion of Bi-CG step
//      y = M\q;
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
//    beta_n_ = beta_d_
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
//    alpha_ = beta_n_ / ( vrs_ );
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

//----------------------------------------------------------------------

EnzoMethodGravityBiCGStab::EnzoMethodGravityBiCGStab 
(const FieldDescr * field_descr,
 int iter_max, double res_tol) 
  : Method(), iter_max_(iter_max), res_tol_(res_tol)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | iter_max_;
  p | res_tol_;

}

//----------------------------------------------------------------------

void EnzoMethodGravityBiCGStab::compute ( Block * block) throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  Field field = enzo_block->data()->field();

  // density field

  const int id = field.field_id("density");
  int gd3[3];
  field.ghost_depth(id,&gd3[0],&gd3[1],&gd3[2]);
  int nd3[3];
  field.size(&nd3[0],&nd3[1],&nd3[2]);
  int md3[3] = {nd3[0] > 1 ? nd3[0]+2*gd3[0] : 1,
		nd3[1] > 1 ? nd3[1]+2*gd3[1] : 1,
		nd3[2] > 1 ? nd3[2]+2*gd3[2] : 1};
  void * density = field.values(id);

  // gravitational potential field
  const int ip = field.field_id("potential");
  int gp3[3];
  field.ghost_depth(ip,&gp3[0],&gp3[1],&gp3[2]);
  int np3[3];
  field.size(&np3[0],&np3[1],&np3[2]);
  int mp3[3] = {np3[0] > 1 ? np3[0]+2*gp3[0] : 1,
		np3[1] > 1 ? np3[1]+2*gp3[1] : 1,
		np3[2] > 1 ? np3[2]+2*gp3[2] : 1};
  void * potential = field.values(ip);

  // precision
  const int p = field.precision(id);

  if      (p == precision_single)
    compute_
      ( block, 
	(float*) density, md3,nd3, 
	(float*) potential, mp3,np3);
  else if (p == precision_double)
    compute_
      ( block, 
	(double*) density, md3,nd3, 
	(double*) potential, mp3,np3);
  else if (p == precision_quadruple)
    compute_
      ( block, 
	(long double*) density, md3,nd3, 
	(long double*) potential, mp3,np3);
  else 
    ERROR1("EnzoMethodGravityBiCGStab()", "precision %d not recognized", p);
}

//======================================================================

template <class T>
void EnzoMethodGravityBiCGStab::compute_ 
(Block * block,
 T * density,   int md3[3], int nd3[3],
 T * potential, int mp3[3], int np3[3]) const throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  printf ("potential size %d %d %d  dim %d %d %d\n",
	  np3[0],np3[1],np3[2],
	  mp3[0],mp3[1],mp3[2]);

  printf ("density size %d %d %d  dim %d %d %d\n",
	  nd3[0],nd3[1],nd3[2],
	  md3[0],md3[1],md3[2]);

  enzo_block->compute_done();

}
