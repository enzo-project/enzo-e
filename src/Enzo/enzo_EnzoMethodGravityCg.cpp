// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityCg class
///
/// function [x] = conjgrad(A,b,x)
///     r=b-A*x;
///     p=r;
///     rsold=r'*r;
///  
///     for i=1:1e6
///         Ap=A*p;
///         alpha=rsold/(p'*Ap);
///         x=x+alpha*p;
///         r=r-alpha*Ap;
///         rsnew=r'*r;
///         if sqrt(rsnew)<1e-10
///               break;
///         end
///         p=r+rsnew/rsold*p;
///         rsold=rsnew;
///     end
/// end
/// nabla ^ 2 (potential) = 4 pi G density
///
///======================================================================
///
/// nabla ^ 2 (potential) = 4 pi G density
///
/// cg_begin:
///
///    B = 4 * PI * G * density
///
///    R = MATVEC (A,X) ==> cg_begin_1
///
/// cg_begin_1:
///
///    R = B - R;
///    P = R
///    iter_ = 0
///    rr_ = DOT(R,R) ==> cg_loop_1
///
/// cg_loop_1:
///
///    if (iter_ > iter_max_) ==> cg_end (return_error_max_iter_reached)
///    AP = MATVEC (A,P) ==> cg_loop_2
///
/// cg_loop_2:
///
///    pap_ = DOT(P, AP) ==> cg_loop_3
///
/// cg_loop_3:
///
///    alpha_ = rr_ / pap_
///    X = X + alpha_ * P;
///    R = R - alpha_ * AP;
///    rr_new_ = DOT(R,R) ==> cg_loop_4
///
/// cg_loop_4:
///
///    if (rr_new_ < resid_tol*resid_tol) ==> cg_end(return_converged)
///
///    P = R + rr_new_ / rr_ * P;
/// 
///    rr_ = rr_new_
///    iter = iter + 1
///    ==> cg_loop_1()
///
/// cg_end (return):
///
///    if (return == return_converged) {
///       potential = X
///       ==> cg_exit()
///    } else {
///       ERROR (return-)
///    }
///
/// cg_exit()
///
///    deallocate
///    NEXT()


#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

//----------------------------------------------------------------------

EnzoMethodGravityCg::EnzoMethodGravityCg (int iter_max, double res_tol) 
  : Method(), 
    iter_max_(iter_max), 
    res_tol_(res_tol),
    precision_(precision_unknown),
    /// Input / output vectors
    enzo_block_(0),
    density_(0),
    potential_(0),
    /// CG temporary vectors
    B_(0),
    X_(0),
    R_(0),
    P_(0),
    W_(0),
    AP_(0),
    /// vector attributes
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    /// CG scalars
    iter_(0),
    alpha_(0),
    pap_(0),
    rr_(0),
    rr_new_(0)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | iter_max_;
  p | res_tol_;

  WARNING("EnzoMethodGravityCg::pup()", "skipping transient attributes");

}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::compute ( CommBlock * comm_block) throw()
{

  enzo_block_ = static_cast<EnzoBlock*> (comm_block);
  Field field = enzo_block_->block()->field();

  const int id = field.field_id("density");
  density_ = field.values(id);
  const int ip = field.field_id("potential");
  void * potential = field.values(ip);

  // ASSUME density AND potential FIELD ATTRIBUTES THE SAME

  field.ghosts    (id,&gx_,&gy_,&gz_);
  field.size         (&nx_,&ny_,&nz_);
  field.dimensions(id,&mx_,&my_,&mz_);

  // store precision
  precision_ = field.precision(id);

  if      (precision_ == precision_single)    compute_<float>();
  else if (precision_ == precision_double)    compute_<double>();
  else if (precision_ == precision_quadruple) compute_<long double>();
  else 
    ERROR1("EnzoMethodGravityCg()", "precision %d not recognized", precision_);
}

//======================================================================

template <class T>
void EnzoMethodGravityCg::compute_ () throw()
/// cg_begin:
///
///    B = 4 * PI * G * density
///
///    R = MATVEC (A,X) ==> cg_begin_1
///
{

  printf ("%s:%d potential size %d %d %d  dim %d %d %d\n",
	  __FILE__,__LINE__,
	  nx_,ny_,nz_,
	  mx_,my_,mz_);

  // local 
  T * density   = (T *) density_;
  T * potential = (T *) potential_;
  T * B         = (T *) B_;
  T * X         = (T *) X_;
  T * R         = (T *) R_;
  T * P         = (T *) P_;
  T * W         = (T *) W_;
  T * AP        = (T *) AP_;

  /// Allocate vectors
  /// ERROR: these should be associated with CommBlock blocks not Simulation Method
  B =  new T[mx_*my_*mz_];
  X =  new T[mx_*my_*mz_];
  R =  new T[mx_*my_*mz_];
  P =  new T[mx_*my_*mz_];
  W =  new T[mx_*my_*mz_];
  AP = new T[mx_*my_*mz_];

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] = 0.0;
	B[i] = 4.0 * (cello::pi) * (cello::G_cgs) * density[i];
	P[i] = R[i] = B[i];
      }
    }
  }

  // rr_ = DOT(R,R)

  // local contribution
  double b_sum = 0.0;
  double r_sum = 0.0;
  const int i0 = gx_ + mx_*(gy_ + my_*gz_);
  double rr = 0.0;
  for (int iz=0; iz<nz_; iz++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + ix + mx_*(iy + my_*iz);
	b_sum += B[i];
	r_sum += R[i];
	rr += R[i]*R[i];
      }
    }
  }
  printf ("%s:%d b_sum = %g\n",__FILE__,__LINE__,b_sum);
  printf ("%s:%d r_sum = %g\n",__FILE__,__LINE__,r_sum);

  iter_ = 0;

  CkCallback callback(CkIndex_EnzoBlock::r_method_gravity_cg_1(NULL),
		      enzo_block_->proxy_array());
  printf ("%s:%d rr local = %g\n",__FILE__,__LINE__,rr);
  enzo_block_->contribute (sizeof(double), &rr, CkReduction::sum_double, callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_cg_1 (CkReductionMsg * msg)
{
  double rr = ((double*)msg->getData())[0];
  delete msg;
  printf ("%s:%d rr global = %g\n",__FILE__,__LINE__,rr);
  EnzoMethodGravityCg * method_cg = static_cast<EnzoMethodGravityCg*> (this->method());
  method_cg->cg_loop_1(rr);
}

void EnzoMethodGravityCg::cg_loop_1 (double rr) throw()
{
  rr_ = rr;
  printf ("%s:%d rr_ global = %g\n",__FILE__,__LINE__,rr_);
  if (precision_ == precision_single)    cg_exit_<float>();
  if (precision_ == precision_double)    cg_exit_<double>();
  if (precision_ == precision_quadruple) cg_exit_<long double>();
}

/// nabla ^ 2 (potential) = 4 pi G density
///
/// cg_begin_1:
///
///    R = B - R;
///    P = R
///    iter_ = 0
///    rr_ = DOT(R,R) ==> cg_loop_1
///
/// cg_loop_1:
///
///    if (iter_ > iter_max_) ==> cg_end (return_error_max_iter_reached)
///    AP = MATVEC (A,P) ==> cg_loop_2
///
/// cg_loop_2:
///
///    pap_ = DOT(P, AP) ==> cg_loop_3
///
/// cg_loop_3:
///
///    alpha_ = rr_ / pap_
///    X = X + alpha_ * P;
///    R = R - alpha_ * AP;
///    rr_new_ = DOT(R,R) ==> cg_loop_4
///
/// cg_loop_4:
///
///    if (rr_new_ < resid_tol*resid_tol) ==> cg_end(return_converged)
///
///    P = R + rr_new_ / rr_ * P;
/// 
///    rr_ = rr_new_
///    iter = iter + 1
///    ==> cg_loop_1()
///
/// cg_end (return):
///    if (return == return_converged) {
///       potential = X
///       ==> cg_exit()
///    } else {
///       ERROR (return-)
///    }
///
template <class T>
void EnzoMethodGravityCg::cg_exit_() throw()
{
  /// deallocate vectors

  delete [] (T *) B_;  B_ = 0;
  delete [] (T *) X_;  X_ = 0;
  delete [] (T *) R_;  R_ = 0;
  delete [] (T *) P_;  P_ = 0;
  delete [] (T *) W_;  W_ = 0;
  delete [] (T *) AP_; AP_ = 0;
  
}

