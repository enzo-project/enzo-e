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

//----------------------------------------------------------------------

EnzoMethodGravityCg::EnzoMethodGravityCg (int iter_max, double res_tol) 
  : Method(), 
    iter_max_(iter_max), 
    res_tol_(res_tol),
    precision_(precision_unknown),
    /// Input / output vectors
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

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);
  Field field = enzo_block->block()->field();

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
{

  printf ("potential size %d %d %d  dim %d %d %d\n",
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
  B =  new T[nx_*ny_*nz_];
  X =  new T[nx_*ny_*nz_];
  R =  new T[nx_*ny_*nz_];
  P =  new T[nx_*ny_*nz_];
  W =  new T[nx_*ny_*nz_];
  AP = new T[nx_*ny_*nz_];

  for (int iz=0; iz<nz_; iz++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = ix + mx_*(iy + my_*iz); // Fields
	int j = ix + nx_*(iy + ny_*iz); // allocated vectors
	B[j] = density[i];
      }
    }
  }
  cg_exit_<T>();
}

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

  delete [] (T *) B_;
  delete [] (T *) X_;
  delete [] (T *) R_;
  delete [] (T *) P_;
  delete [] (T *) W_;
  delete [] (T *) AP_;
  
}

