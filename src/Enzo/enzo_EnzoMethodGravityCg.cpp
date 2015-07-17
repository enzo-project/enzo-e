// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityCg class
///
/// Preconditioned Conjugate Gradient Method
///
/// - A*X = B, preconditioner M
/// - nabla ^ 2 (potential) = 4 pi G density
///
///     - X = initial guess
///     - B = right-hand side
///     - R = B - A*X
///     - solve(M*Z = R)
///     - D = Z
///     - shift (B)
///
///     - do
///        -  Y = A*D;
///        -  rz = dot(R,Z)
///        -  dy = dot(D,Y)
///
///        -  a = rz / dy;
///        -  X = X + a*D;
///        -  R = R - a*Y;
///        -  solve(M*Z = R)
///
///        -  rz2 = dot(R,Z)
///        -  b = rz2 / rz;
///        -  D = Z + b*D;
///        -  rz = rz2;
///     - while (! converged())
///
/// Required Fields
/// 
/// - B                          linear system right-hand side
/// - D                          CG search direction vector: unpreconditioned
/// - R                          residual vector of unpreconditioned linear system
/// - X                          linear system right-hand solution
/// - Y                          Krylov subspace vector
/// - Z                          residual of preconditioned linear system
/// - potential                  computed gravitational potential
/// - density                    density field
/// - acceleration_x             acceleration along X-axis
/// - acceleration_y [rank >= 2] acceleration along Y-axis
/// - acceleration_z [rank >= 3] acceleration along Z-axis


#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

/* #define NEW_REFRESH */
/* #define DEBUG_GRAVITY */

#ifndef NEW_REFRESH
#  define OLD_REFRESH
#endif


#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#ifdef DEBUG_METHOD
#   define TRACE_METHOD					\
  CkPrintf ("%s:%d TRACE_METHOD\n",__FILE__,__LINE__);	\
  fflush(stdout);
#else
#   define TRACE_METHOD /*  */ 
#endif
//----------------------------------------------------------------------

EnzoMethodGravityCg::EnzoMethodGravityCg 
(const FieldDescr * field_descr, int rank,
 double grav_const, int iter_max, double res_tol, int monitor_iter,
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
    monitor_iter_(monitor_iter),
    rr0_(0),
    rr_min_(0),rr_max_(0),
    idensity_(0),  ipotential_(0),
    ib_(0), ix_(0), ir_(0), id_(0), iy_(0), iz_(0),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    rr_(0.0), rz_(0.0), rz2_(0.0), dy_(0.0), bs_(0.0), bc_(0.0)
#ifdef OLD_REFRESH
  , id_refresh_matvec_(-1)
#endif
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

  /// Initialize default Refresh

  const int num_fields = field_descr->field_count();

  field_descr->ghost_depth    (idensity_,&gx_,&gy_,&gz_);

  const int ir = add_refresh(1,rank-1,neighbor_leaf,sync_barrier);
  //  refresh(ir)->add_field(idensity_);
  refresh(ir)->add_all_fields(num_fields);

  /// Initialize matvec Refresh

#ifdef OLD_REFRESH
  id_refresh_matvec_ = add_refresh(1,rank-1,neighbor_leaf,sync_barrier);
  refresh(id_refresh_matvec_)->add_all_fields(num_fields);
  //  refresh(id_refresh_matvec_)->add_field(ir_);
#endif

}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::compute ( Block * block) throw()
{
  
  Field field = block->data()->field();

  field.size                (&nx_,&ny_,&nz_);
  field.dimensions(idensity_,&mx_,&my_,&mz_);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

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

  iter_ = 0;

  Data * data = enzo_block->data();
  Field field = data->field();

  if (enzo_block->is_leaf()) {

    ///   - X = 0
    ///   - B = -h^2 * 4 * PI * G * density
    ///   - R = P = B ( residual with X = 0);

    T * density = (T*) field.values(idensity_);
    
    T * B = (T*) field.values(ib_);
    T * X = (T*) field.values(ix_);
    T * R = (T*) field.values(ir_);

    const int ix0 = 0;
    const int iy0 = 0;
    const int iz0 = 0;

    for (int iz=iz0; iz<mz_-iz0; iz++) {
      for (int iy=iy0; iy<my_-iy0; iy++) {
	for (int ix=ix0; ix<mx_-ix0; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  X[i] = 0.0;
	  B[i] = - 4.0 * (cello::pi) * grav_const_ * density[i];
	  R[i] = B[i];
	}
      }
    }

    M_->matvec(id_,ir_,enzo_block);
    M_->matvec(iz_,ir_,enzo_block);
  }

  long double reduce[3];

  if (enzo_block->is_leaf()) {

    T * B = (T*) field.values(ib_);
    T * R = (T*) field.values(ir_);
    reduce[0] = dot_(R,R);
    reduce[1] = sum_(B);
    reduce[2] = count_();
    
  } else {

    reduce[0] = 0.0;
    reduce[1] = 0.0;
    reduce[2] = 0.0;

  }

  CkCallback callback(CkIndex_EnzoBlock::r_cg_loop_0a<T>(NULL), 
		      enzo_block->proxy_array());

#ifdef DEBUG_GRAVITY
  printf ("%s:%d %s DEBUG_GRAVITY calling contribute\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif
	  
  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  r_method_gravity_cg_type, 
			  callback);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_cg_loop_0a (CkReductionMsg * msg)
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> refresh P for AP = MATVEC (A,P)
{
  TRACE_METHOD;

  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  long double * data = (long double *) msg->getData();

  method->set_rr( data[0] );
  method->set_bs( data[1] );
  method->set_bc( data[2] );

  delete msg;
  
  // Refresh if Block is a leaf

#ifdef OLD_REFRESH
  method->refresh(1)->set_active(is_leaf());
  refresh_enter(CkIndex_EnzoBlock::r_enzo_matvec(NULL),
		method->refresh(1));
#endif
#ifdef NEW_REFRESH
  Refresh refresh (4,0,neighbor_level, sync_face);
  refresh.set_active(is_leaf());
  refresh.add_all_fields(this->data()->field().field_count());
  refresh_enter(CkIndex_EnzoBlock::r_enzo_matvec(NULL),
		&refresh);
#endif
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_cg_loop_0b (CkReductionMsg * msg)
/// ==> refresh P for AP = MATVEC (A,P)
{

  EnzoMethodGravityCg * method = 
    static_cast<EnzoMethodGravityCg*> (this->method());

  method->set_iter ( ((int*)msg->getData())[0] );

  // Refresh if Block is a leaf, then continue with enzo_matvec_()

#ifdef OLD_REFRESH
  method->refresh(1)->set_active(is_leaf());
  refresh_enter(CkIndex_EnzoBlock::r_enzo_matvec(NULL),
		method->refresh(1));
#endif

#ifdef NEW_REFRESH
  Refresh refresh (4,0,neighbor_level, sync_face);
  refresh.set_active(is_leaf());
  refresh.add_all_fields(this->data()->field().field_count());
  refresh_enter(CkIndex_EnzoBlock::r_enzo_matvec(NULL),
		&refresh);
#endif
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

  Data * data = enzo_block->data();
  Field field = data->field();

  if (enzo_block->is_leaf()) {

    cello::check(rr_,"rr_",__FILE__,__LINE__);
    cello::check(bs_,"bs_",__FILE__,__LINE__);
    cello::check(bc_,"bc_",__FILE__,__LINE__);

    T * B  = (T*) field.values(ib_);
    T * R  = (T*) field.values(ir_);

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
  }

  long double reduce[3] = {0.0, 0.0, 0.0};

  if (enzo_block->is_leaf()) {

    T * R  = (T*) field.values(ir_);
    reduce[0] = dot_(R,R);

  } else {

    reduce[0] = 0.0;

  } 

  CkCallback callback(CkIndex_EnzoBlock::r_cg_shift_1<T>(NULL), 
		enzo_block->proxy_array());

#ifdef DEBUG_GRAVITY
  printf ("%s:%d %s DEBUG_GRAVITY calling contribute\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif

  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  r_method_gravity_cg_type, 
			  callback);
    
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

    Monitor * monitor = enzo_block->simulation()->monitor();

    if (iter_ == 0) {
      monitor->print ("Enzo", "CG iter %04d  rr0 %g",
		      iter_,(double)(rr0_));
    }

    if (monitor_iter_ && (iter_ % monitor_iter_) == 0 ) {
      monitor_output_ (enzo_block);
    }
  }

  if (rr_ / rr0_ < res_tol_) {

    cg_end<T> (enzo_block,return_converged);

  } else if (iter_ >= iter_max_)  {

    cg_end<T>(enzo_block,return_error_max_iter_reached);

  } else {
    
    Data * data = enzo_block->data();
    Field field = data->field();

    if (enzo_block->is_leaf()) {

      double hx,hy,hz;
      data->field_cell_width(&hx,&hy,&hz);

      A_->matvec(iy_,id_,enzo_block);

    }

    long double reduce[3];

    if (enzo_block->is_leaf()) {

      T * D = (T*) field.values(id_);
      T * Y = (T*) field.values(iy_);
      T * R = (T*) field.values(ir_);
      T * Z = (T*) field.values(iz_);

      reduce[0] = dot_(R,R);
      reduce[1] = dot_(R,Z);
      reduce[2] = dot_(D,Y);

    } else {

      reduce[0] = 0.0;
      reduce[1] = 0.0;
      reduce[2] = 0.0;
    }

    CkCallback callback(CkIndex_EnzoBlock::r_cg_loop_3<T>(NULL), 
		  enzo_block->proxy_array());

#ifdef DEBUG_GRAVITY
  printf ("%s:%d %s DEBUG_GRAVITY calling contribute\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif

    enzo_block->contribute (3*sizeof(long double), &reduce, 
			    r_method_gravity_cg_type,
			    callback);
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

  cello::check(rr_,"rr_",__FILE__,__LINE__);
  cello::check(rz_,"rz_",__FILE__,__LINE__);
  cello::check(dy_,"dy_",__FILE__,__LINE__);

  Data * data = enzo_block->data();
  Field field = data->field();

  if (enzo_block->is_leaf()) {

    T * X = (T*) field.values(ix_);
    T * D = (T*) field.values(id_);
    T * R = (T*) field.values(ir_);
    T * Y = (T*) field.values(iy_);

    T a = rz_ / dy_;

    cello::check(a,"a",__FILE__,__LINE__);

    zaxpy_ (X,  a ,D,X);
    zaxpy_ (R, -a, Y,R);

    M_->matvec(iz_,ir_,enzo_block);

  }

  long double reduce[3];

  if (enzo_block->is_leaf()) {

    T * X = (T*) field.values(ix_);
    T * R = (T*) field.values(ir_);
    T * Z = (T*) field.values(iz_);
    reduce[0] = dot_(R,Z);
    reduce[1] = sum_(R);
    reduce[2] = sum_(X);

  } else {

    reduce[0] = 0.0;
    reduce[1] = 0.0;
    reduce[2] = 0.0;

  }

  CkCallback callback(CkIndex_EnzoBlock::r_cg_loop_5<T>(NULL), 
		      enzo_block->proxy_array());

#ifdef DEBUG_GRAVITY
  printf ("%s:%d %s DEBUG_GRAVITY calling contribute\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif

  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  r_method_gravity_cg_type, 
			  callback);
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

  cello::check(rz2_,"rz2_",__FILE__,__LINE__);
  cello::check(rs_,"rs_",__FILE__,__LINE__);
  cello::check(xs_,"xs_",__FILE__,__LINE__);

  Field field = enzo_block->data()->field();

  if (enzo_block->is_leaf()) {

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

    cello::check(b,"b",__FILE__,__LINE__);

    zaxpy_ (D,b,D,Z);
  }

  int iter = iter_ + 1;

  CkCallback callback(CkIndex_EnzoBlock::r_cg_loop_0b<T>(NULL), 
		enzo_block->proxy_array());
    
  enzo_block->contribute (sizeof(int), &iter, 
			  CkReduction::max_int, callback);

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
  if (enzo_block->is_leaf()) {

    Data * data = enzo_block->data();
    Field field = data->field();

    T * X         = (T*) field.values(ix_);
    T * potential = (T*) field.values(ipotential_);

    if (enzo_block->index().is_root()) {
      monitor_output_ (enzo_block);
    }

    copy_(potential,X,mx_,my_,mz_);

    bool symmetric;
    int order;
    EnzoComputeAcceleration compute_acceleration (field.field_descr(),
						  rank_, symmetric = true,
						  order=2);
    compute_acceleration.compute(enzo_block);

  }

  enzo_block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::monitor_output_(EnzoBlock * enzo_block) throw()
{

  Monitor * monitor = enzo_block->simulation()->monitor();

  monitor->print ("Enzo", "CG iter %04d  rr %g [%g %g]",
		  iter_,
		  (double)(rr_    / rr0_),
		  (double)(rr_min_/ rr0_),
		  (double)(rr_max_/ rr0_));
}

//----------------------------------------------------------------------

void EnzoMethodGravityCg::cg_exit_() throw()
/// deallocate temporary vectors
{
}

//======================================================================

template <class T>
long double EnzoMethodGravityCg::dot_ (const T * X, const T * Y) const throw()
{
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

int EnzoMethodGravityCg::count_ () const throw()
{
  return nx_*ny_*nz_;
}

//----------------------------------------------------------------------
