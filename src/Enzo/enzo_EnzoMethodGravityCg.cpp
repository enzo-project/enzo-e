// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityCg class
//
// Preconditioned Conjugate Gradient Method
//
// A*X = B, preconditioner M
// nabla ^ 2 (potential) = 4 pi G density
//
//     X = initial guess
//     B = right-hand side
//     R = B - A*X
//     solve(M*Z = R)
//     D = Z
//     shift (B)

//     do
//         Y = A*D;
//         rz = dot(R,Z)
//         dy = dot(D,Y)

//         a = rz / dy;
//         X = X + a*D;
//         R = R - a*Y;
//         solve(M*Z = R)
//
//         rz2 = dot(R,Z)
//         b = rz2 / rz;
//         D = Z + b*D;
//         rz = rz2;
//     while (! converged())
//

#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#define CHECK(VALUE,STRING)			\
  if (std::fpclassify(VALUE) == FP_ZERO) {	\
    printf ("ERROR: %s:%d %s zero\n",		\
	    __FILE__,__LINE__,STRING);		\
  }						\
  if (std::fpclassify(VALUE) == FP_NAN) {	\
    printf ("ERROR: %s:%d %s nan\n",		\
	    __FILE__,__LINE__,STRING);		\
  }						\
  if (std::fpclassify(VALUE) == FP_INFINITE) {	\
    printf ("ERROR: %s:%d %s inf\n",		\
	    __FILE__,__LINE__,STRING);		\
  }						\

//----------------------------------------------------------------------

EnzoMethodGravityCg::EnzoMethodGravityCg 
(FieldDescr * field_descr, int rank,
 double grav_const, int iter_max, double res_tol,
 bool is_singular,
 bool diag_precon) 
  : Method(), 
    A_(new EnzoMatrixLaplace),
    M_(NULL),
    is_singular_(is_singular),
    rank_(rank),
    grav_const_(grav_const),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    rr0_(0),
    rr_min_(0),rr_max_(0),
    idensity_(0),  ipotential_(0),
    ib_(0), ix_(0), ir_(0), id_(0), iy_(0), iz_(0),
    is_leaf_(false),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    rr_(0.0), rz_(0.0), rz2_(0.0), dy_(0.0), bs_(0.0), bc_(0.0)
{

  M_ = (diag_precon) ? (Matrix *)(new EnzoMatrixDiagonal) 
    :                  (Matrix *)(new EnzoMatrixIdentity);

  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  id_ = field_descr->field_id("D");
  ir_ = field_descr->field_id("R");
  ix_ = field_descr->field_id("X");
  iy_ = field_descr->field_id("Y");
  iz_ = field_descr->field_id("Z");

}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::compute ( Block * block) throw()
{
  
  set_leaf(block);

  Field field = block->data()->field();

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  field.size                (&nx_,&ny_,&nz_);
  field.dimensions(idensity_,&mx_,&my_,&mz_);
  field.ghosts    (idensity_,&gx_,&gy_,&gz_);

  int precision = field.precision(idensity_);

  if      (precision == precision_single)    compute_<float>      (enzo_block);
  else if (precision == precision_double)    compute_<double>     (enzo_block);
  else if (precision == precision_quadruple) compute_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityCg()", "precision %d not recognized", precision);
}

//======================================================================

extern CkReduction::reducerType r_method_gravity_cg_type;

template <class T>
void EnzoMethodGravityCg::compute_ (EnzoBlock * enzo_block) throw()
//     X = initial guess
//     B = right-hand side
//     R = B - A*X
//     solve(M*Z = R)
//     D = Z
//     shift (B)
{

  set_leaf(enzo_block);

  iter_ = 0;

  Data * data = enzo_block->data();
  Field field = data->field();

  T * density = (T*) field.values(idensity_);
    
  T * B = (T*) field.values(ib_);
  T * X = (T*) field.values(ix_);
  T * R = (T*) field.values(ir_);
  T * D = (T*) field.values(id_);
  T * Y = (T*) field.values(iy_);
  T * Z = (T*) field.values(iz_);

  if (is_leaf_) {

    ///   - X = 0
    ///   - B = -h^2 * 4 * PI * G * density
    ///   - R = P = B ( residual with X = 0);

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

  CkCallback cb(CkIndex_EnzoBlock::r_cg_loop_0a<T>(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  r_method_gravity_cg_type, 
			  cb);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_cg_loop_0a (CkReductionMsg * msg)
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> refresh P for AP = MATVEC (A,P)
{
  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  long double * data = (long double *) msg->getData();

  method->set_rr( data[0] );
  method->set_bs( data[1] );
  method->set_bc( data[2] );

  delete msg;

  refresh_sync_  = "contribute";
  refresh_phase_ = phase_enzo_matvec;

  index_refresh_ = method->index_refresh(1);

  control_sync(phase_refresh_enter,"neighbor");
  
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_cg_loop_0b (CkReductionMsg * msg)
/// ==> refresh P for AP = MATVEC (A,P)
{

  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  method->set_iter ( ((int*)msg->getData())[0] );

  delete msg;

  refresh_sync_  = "contribute";
  refresh_phase_ = phase_enzo_matvec;

  index_refresh_ = method->index_refresh(1);

  control_sync(phase_refresh_enter,"neighbor");
}

//----------------------------------------------------------------------

void EnzoBlock::enzo_matvec_()
{
  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (this);

  if      (precision == precision_single)    
    method->cg_shift_1<float>(enzo_block);
  else if (precision == precision_double)    
    method->cg_shift_1<double>(enzo_block);
  else if (precision == precision_quadruple) 
    method->cg_shift_1<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityCg()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityCg::cg_shift_1 (EnzoBlock * enzo_block) throw()
{

  CHECK(rr_,"rr_");
  CHECK(bs_,"bs_");
  CHECK(bc_,"bc_");

  set_leaf(enzo_block);

  Data * data = enzo_block->data();
  Field field = data->field();

  T * B  = (T*) field.values(ib_);
  T * R  = (T*) field.values(ir_);
  T * D  = (T*) field.values(id_);
  T * Z  = (T*) field.values(iz_);

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

  CkCallback cb(CkIndex_EnzoBlock::r_cg_shift_1<T>(NULL), 
		enzo_block->proxy_array());

  enzo_block->contribute (sizeof(long double), &reduce, 
			  r_method_gravity_cg_type, 
			  cb);
    
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_cg_shift_1 (CkReductionMsg * msg)
{
  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  method->set_rr( ((long double*)msg->getData())[0] );

  delete msg;

  method -> cg_loop_2<T>(this);

}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityCg::cg_loop_2 (EnzoBlock * enzo_block) throw()
{
  set_leaf(enzo_block);

  CHECK(rr_,"rr_");

  if (iter_ == 0) {
    rr0_ = rr_;
    rr_min_ = rr_;
    rr_max_ = rr_;
  } else {
    rr_min_ = std::min(rr_min_,rr_);
    rr_max_ = std::max(rr_max_,rr_);
  }

  if (enzo_block->index().is_root()) {
    printf ("%s:%d cg_end iter %d rr0 %g rr_min %g rr_max %g rr %g\n",
	    __FILE__,__LINE__,iter_,
	    (double)rr0_,(double)rr_min_,(double)rr_max_,(double)rr_);
  }

  if (rr_ / rr0_ < res_tol_) {

    cg_end<T> (enzo_block,return_converged);

  } else if (iter_ >= iter_max_)  {

    cg_end<T>(enzo_block,return_error_max_iter_reached);

  } else {

    Data * data = enzo_block->data();
    Field field = data->field();

    T * D = (T*) field.values(id_);
    T * Y = (T*) field.values(iy_);
    T * R = (T*) field.values(ir_);
    T * Z = (T*) field.values(iz_);

    double hx,hy,hz;
    data->field_cell_width(&hx,&hy,&hz);

    A_->matvec(iy_,id_,enzo_block);

    long double reduce[3];

    reduce[0] = dot_(R,R);
    reduce[1] = dot_(R,Z);
    reduce[2] = dot_(D,Y);

    CkCallback cb(CkIndex_EnzoBlock::r_cg_loop_3<T>(NULL), 
		  enzo_block->proxy_array());

    enzo_block->contribute (3*sizeof(long double), &reduce, 
			    r_method_gravity_cg_type,
			    cb);
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_cg_loop_3 (CkReductionMsg * msg)
{
  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  long double * data = (long double *) msg->getData();

  method->set_rr(data[0]);
  method->set_rz(data[1]);
  method->set_dy(data[2]);
  
  delete msg;

  method -> cg_loop_4<T>(this);

}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityCg::cg_loop_4 (EnzoBlock * enzo_block) throw ()
//  a = rz / dy;
//  X = X + a*D;
//  R = R - a*Y;
//  solve(M*Z = R)
//  rz2 = dot(R,Z)
//  b = rz2 / rz;
//  D = Z + b*D;
//  rz = rz2;
{

  set_leaf(enzo_block);

  CHECK(rr_,"rr_");
  CHECK(rz_,"rz_");
  CHECK(dy_,"dy_");

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X = (T*) field.values(ix_);
  T * D = (T*) field.values(id_);
  T * R = (T*) field.values(ir_);
  T * Y = (T*) field.values(iy_);
  T * Z = (T*) field.values(iz_);

  T a = rz_ / dy_;

  CHECK(a,"a");

  zaxpy_ (X,  a ,D,X);
  zaxpy_ (R, -a, Y,R);

  M_->matvec(iz_,ir_,enzo_block);

  long double reduce[3];

  reduce[0] = dot_(R,Z);
  reduce[1] = sum_(R);
  reduce[2] = sum_(X);

  CkCallback cb(CkIndex_EnzoBlock::r_cg_loop_5<T>(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  r_method_gravity_cg_type, 
			  cb);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_cg_loop_5 (CkReductionMsg * msg)
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> cg_loop_6
{
  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  long double * data = (long double *) msg->getData();

  method->set_rz2(data[0]);
  method->set_rs (data[1]);
  method->set_xs (data[2]);

  delete msg;

  method -> cg_loop_6<T>(this);

}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityCg::cg_loop_6 (EnzoBlock * enzo_block) throw ()
//  rz2 = dot(R,Z)
//  b = rz2 / rz;
//  D = Z + b*D;
//  rz = rz2;
{

  set_leaf(enzo_block);

  CHECK(rz2_,"rz2_");
  CHECK(rs_,"rs_");
  CHECK(xs_,"xs_");

  Field field = enzo_block->data()->field();

  if (is_singular_)  {

    // shift rhs B by projection of B onto e: B~ <== B - (e*eT)/(eT*e) b
    // eT*e == n === zone count (bc)
    // eT*b == sum_i=1,n B[i]

    T * X  = (T*) field.values(ix_);
    T * R  = (T*) field.values(ir_);
    shift_ (X,T(-xs_/bc_),X);
    shift_ (R,T(-rs_/bc_),R);
  }

  T * D  = (T*) field.values(id_);
  T * Z  = (T*) field.values(iz_);

  T b = rz2_ / rz_;

  CHECK(b,"b");

  zaxpy_ (D,b,D,Z);

  int iter = iter_ + 1;

  CkCallback cb(CkIndex_EnzoBlock::r_cg_loop_0b<T>(NULL), 
		enzo_block->proxy_array());
    
  enzo_block->contribute (sizeof(int), &iter, 
			  CkReduction::max_int, cb);

}

//----------------------------------------------------------------------


template <class T>
void EnzoMethodGravityCg::cg_end (EnzoBlock * enzo_block,int retval) throw ()
///    if (return == return_converged) {
///       potential = X
///       ==> cg_exit()
///    } else {
///       ERROR (return-)
///    }
{
  set_leaf(enzo_block);

  enzo_block->clear_refresh();

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X         = (T*) field.values(ix_);
  T * potential = (T*) field.values(ipotential_);

  if (enzo_block->index().is_root()) {
    printf ("%s:%d cg_end iter %d rr0 %g rr_min %g rr_max %g rr %g\n",
	    __FILE__,__LINE__,iter_,
	    (double)rr0_,(double)rr_min_,(double)rr_max_,(double)rr_);
  }

  copy_(potential,X);

  bool symmetric;
  int order;
  EnzoComputeAcceleration compute_acceleration (field.field_descr(),
						rank_, symmetric = true,
						order=2);
  compute_acceleration.compute(enzo_block);
  enzo_block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::cg_exit_() throw()
/// deallocate temporary vectors
{
}

//======================================================================

template <class T>
void EnzoMethodGravityCg::copy_ (T * X, const T * Y) const throw()
{
  if (! is_leaf_ ) return;

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] = Y[i];
      }
    }
  }
}

//======================================================================

template <class T>
long double EnzoMethodGravityCg::dot_ (const T * X, const T * Y) const throw()
{
  if (! is_leaf_ ) return 0.0;

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

template <class T>
void EnzoMethodGravityCg::zaxpy_ (T * Z, double a, const T * X, const T * Y) const throw()
{
  if (! is_leaf_ ) return;

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

template <class T>
long double EnzoMethodGravityCg::sum_ (const T * X) const throw()
{
  if (! is_leaf_ ) return 0.0;

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

template <class T>
void EnzoMethodGravityCg::shift_ (T * X, const T a, const T * Y) const throw()
{
  if (! is_leaf_ ) return;

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

template <class T>
void EnzoMethodGravityCg::scale_ (T * Y, T a, const T * X) const throw()
{
  if (! is_leaf_ ) return;

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

template <class T>
T EnzoMethodGravityCg::count_ (T * X) const throw()
{
  if (! is_leaf_ ) return 0.0;

  T count = 0.0;

  const int i0 = gx_ + mx_*(gy_ + my_*gz_);

  for (int iz=0; iz<nz_; iz++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + (ix + mx_*(iy + my_*iz));
	++count;
      }
    }
  }
  return count;
}

//----------------------------------------------------------------------
