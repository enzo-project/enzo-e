// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityBiCGStab.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-23 16:19:06
/// @brief    Implements the EnzoMethodGravityBiCGStab class
//
// The following is from http://www.netlib.org/templates/matlab/bicgstab.m:
//
// % function [x, error, iter, flag] = bicgstab(A, x, b, M, max_it, tol)
// %
// %  -- Iterative template routine --
// %     Univ. of Tennessee and Oak Ridge National Laboratory
// %     October 1, 1993
// %     Details of this algorithm are described in "Templates for the
// %     Solution of Linear Systems: Building Blocks for Iterative
// %     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
// %     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
// %     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
// %
// %  [x, error, iter, flag] = bicgstab(A, x, b, M, max_it, tol)
// %
// % bicgstab.m solves the linear system Ax=b using the 
// % BiConjugate Gradient Stabilized Method with preconditioning.
// %
// % input   A        REAL matrix
// %         x        REAL initial guess vector
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
//   iter = 0;                                          % initialization
//   flag = 0;
// 
//   bnrm2 = norm( b );
//   if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
// 
//   r = b - A*x;
//   error = norm( r ) / bnrm2;
//   if ( error < tol ) return, end
// 
//   omega  = 1.0;
//   r_tld = r;
// 
//   for iter = 1:max_it,                              % begin iteration
// 
//      rho   = ( r_tld'*r );                          % direction vector
//      if ( rho == 0.0 ) break, end
// 
//      if ( iter > 1 ),
//         beta  = ( rho/rho_1 )*( alpha/omega );
//         p = r + beta*( p - omega*v );
//      else
//         p = r;
//      end
//  
//      p_hat = M \ p;
//      v = A*p_hat;
//      alpha = rho / ( r_tld'*v );
//      s = r - alpha*v;
//      if ( norm(s) < tol ),                          % early convergence check
//         x = x + alpha*p_hat;
//         resid = norm( s ) / bnrm2;
//         break;
//      end
// 
//      s_hat = M \ s;                                 % stabilizer
//      t = A*s_hat;
//      omega = ( t'*s) / ( t'*t );
// 
//      x = x + alpha*p_hat + omega*s_hat;             % update approximation
// 
//      r = s - omega*t;
//      error = norm( r ) / bnrm2;                     % check convergence
//      if ( error <= tol ), break, end
// 
//      if ( omega == 0.0 ), break, end
//      rho_1 = rho;
// 
//   end
// 
//   if ( error <= tol | s <= tol ),                   % converged
//      if ( s <= tol ),
//         error = norm(s) / bnrm2;
//      end
//      flag =  0;
//   elseif ( omega == 0.0 ),                          % breakdown
//      flag = -2;
//   elseif ( rho == 0.0 ),
//      flag = -1;
//   else                                              % no convergence
//      flag = 1;
//   end
// 
// % END bicgstab.m
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
//
//    bnorm_ = NORM( B ) ==> bicgstab_start_1
//
// --------------------
// bicgstab_start_1()
// --------------------
//
//   if ( bnorm_ == 0.0 ) bnorm_ = 1.0;
// 
//   R = MATVEC(A,x) ==> bicgstab_start_2
//
// --------------------
// bicgstab_start_2()
// --------------------
//
//   R = B - R;
//   rnorm_ = NORM( R ) ==> bicgstab_start_3
//
// --------------------
// bicgstab_start_3()
// --------------------
//
//   error = rnorm_ / bnorm_;
//   if ( error < tol ) {
//      ==> bicgstab_loop_end(return_converged_)
//   }
// 
//   R_hat = R;
// 
//   iter = 1
// 
//   ==> bicgstab_loop_begin()
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
//    rho_ = DOT(R_hat, R) ==> bicgstab_loop_1
// 
// --------------------
// bicgstab_loop_1()
// --------------------
//  
//    if ( rho_ == 0.0 ) {
//       bicgstab_loop_end(return_error_rho_eq_0)
//    }
// 
//    if ( iter_ > 1 ) {
//       beta  = ( rho_/rho_prev_ )*( alpha_/omega_ );
//       P = R + beta*( P - omega_*v );
//    } else {
//       P = R;
//    }
//  
//    Y = SOLVE (M, P) ==> bicgstab_loop_2
// 
// --------------------
// bicgstab_loop_2()
// --------------------
//
//    V = MATVEC (A,Y) ==> bicgstab_loop_3
//
// --------------------
// bicgstab_loop_3()
// --------------------
//
//    V = A * Y;
//    rhdv_ = DOT (R_hat,V) ==> bicgstab_loop_4
//
// --------------------
// bicgstab_loop_4()
// --------------------

//    alpha_ = rho / ( rhdv_ );
//    S = R - alpha_*V;
//    Z = SOLVE(M , S) ==> bicgstab_loop_5
//
// --------------------
// bicgstab_loop_5()
// --------------------
//
//    T = MATVEC(A,Z) ==> bicgstab_loop_6
// 
// --------------------
// bicgstab_loop_6()
// --------------------
//    tds_ = DOT(T,S), tdt_ = DOT(T,T) ==> bicgstab_loop_7
//
// --------------------
// bicgstab_loop_7()
// --------------------
//
//    omega_ = ( tds_) / ( tdt_ );
//    if ( omega_ == 0.0 ) {
//       bicgstab_loop_end(return_error_omega_eq_0);
//    }
//    X = X + alpha_*Y + omega_*Z;             % update approximation
//    R = S - omega_*T;
//    rnorm_ = NORM(R) ==> bicgstab_loop_8()
//
// --------------------
// bicgstab_loop_7()
// --------------------
//
//    error = rnorm_ / bnorm_;                     % check convergence
//    if ( error <= tol ) {
//       bicgstab_loop_end(return_converged_);
//    }
// 
//    rho_prev = rho;
//
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
  field.ghosts(id,&gd3[0],&gd3[1],&gd3[2]);
  int nd3[3];
  field.size(&nd3[0],&nd3[1],&nd3[2]);
  int md3[3] = {nd3[0] > 1 ? nd3[0]+2*gd3[0] : 1,
		nd3[1] > 1 ? nd3[1]+2*gd3[1] : 1,
		nd3[2] > 1 ? nd3[2]+2*gd3[2] : 1};
  void * density = field.values(id);

  // gravitational potential field
  const int ip = field.field_id("potential");
  int gp3[3];
  field.ghosts(ip,&gp3[0],&gp3[1],&gp3[2]);
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

  enzo_block->compute_stop();

}
