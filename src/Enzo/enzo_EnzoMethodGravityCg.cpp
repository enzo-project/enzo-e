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

EnzoMethodGravityCg::EnzoMethodGravityCg 
(const FieldDescr * field_descr,
 int iter_max, double res_tol) 
  : Method(), 
    iter_max_(iter_max), 
    res_tol_(res_tol),
    precision_(precision_unknown),
    /// Input / output vectors
    idensity_(0),
    ipotential_(0),
    /// CG temporary vectors
    ib_(0),
    ix_(0),
    ir_(0),
    ip_(0),
    iw_(0),
    iap_(0),
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
  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");
  ib_  = field_descr->field_id("B");
  ix_  = field_descr->field_id("X");
  ir_  = field_descr->field_id("R");
  ip_  = field_descr->field_id("P");
  iw_  = field_descr->field_id("W");
  iap_ = field_descr->field_id("AP");

  WARNING ("EnzoMethodGravityCg::EnzoMethodGravityCg()",
	   "Assuming same ghost depth and precision for all fields");

  field_descr->ghosts    (idensity_,&gx_,&gy_,&gz_);
  precision_ = field_descr->precision(idensity_);

}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | iter_max_;
  p | res_tol_;
  p | precision_;

  p | idensity_;
  p | ipotential_;
  p | ib_;
  p | ix_;
  p | ir_;
  p | ip_;
  p | iw_;
  p | iap_;

  p | nx_;
  p | ny_;
  p | nz_;

  p | mx_;
  p | my_;
  p | mz_;

  p | gx_;
  p | gy_;
  p | gz_;

  WARNING("EnzoMethodGravityCg::pup()", 
	  "skipping transient attributes: iter_,alpha_,pap_,rr_,rr_new_");

}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::compute ( CommBlock * comm_block) throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);
  Field field = enzo_block->block()->field();

  WARNING ("EnzoMethodGravityCg::EnzoMethodGravityCg()",
	   "Assuming same ghost depth and dimensions for all fields");
  field.size                (&nx_,&ny_,&nz_);
  field.dimensions(idensity_,&mx_,&my_,&mz_);

  if      (precision_ == precision_single)    compute_<float>(enzo_block);
  else if (precision_ == precision_double)    compute_<double>(enzo_block);
  else if (precision_ == precision_quadruple) compute_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityCg()", "precision %d not recognized", precision_);
}

//======================================================================

template <class T>
void EnzoMethodGravityCg::compute_ (EnzoBlock * enzo_block) throw()
/// cg_begin:
///
///    B = 4 * PI * G * density
///
///    R = MATVEC (A,X) ==> cg_begin_1
///
{

  Field field = enzo_block->block()->field();

  T * density = (T*) field.values(idensity_);
  T * potential = (T*) field.values(ipotential_);

  T * B = (T*) field.values(ib_);
  T * X = (T*) field.values(ix_);
  T * R = (T*) field.values(ir_);
  T * P = (T*) field.values(ip_);
  T * W = (T*) field.values(iw_);
  T * AP = (T*) field.values(iap_);

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
  const int i0 = gx_ + mx_*(gy_ + my_*gz_);
  double rr = 0.0;
  for (int iz=0; iz<nz_; iz++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + ix + mx_*(iy + my_*iz);
	rr += R[i]*R[i];
      }
    }
  }

  iter_ = 0;

  CkCallback callback(CkIndex_EnzoBlock::r_method_gravity_cg_1(NULL),
		      enzo_block->proxy_array());
  enzo_block->contribute (sizeof(double), &rr, CkReduction::sum_double, callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_cg_1 (CkReductionMsg * msg)
{
  double rr = ((double*)msg->getData())[0];
  delete msg;
  EnzoMethodGravityCg * method_cg = static_cast<EnzoMethodGravityCg*> (this->method());
  method_cg->cg_loop_1(rr);
}

void EnzoMethodGravityCg::cg_loop_1 (double rr) throw()
{
  rr_ = rr;
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

}

