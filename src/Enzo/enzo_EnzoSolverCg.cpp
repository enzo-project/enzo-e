// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-1-08
/// @brief    Implements the CG Krylov iterative linear solver

#include "enzo.hpp"

#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

// #define DEBUG_SOLVER
#ifdef DEBUG_SOLVER
#   define TRACE_SOLVER(solver)						\
  CkPrintf ("%d %s:%d TRACE %s %p\n",CkMyPe(),__FILE__,__LINE__,solver,this); \
  fflush(stdout);
#else
#   define TRACE_SOLVER(solver) /*  */ 
#endif

//----------------------------------------------------------------------

EnzoSolverCg::EnzoSolverCg 
(const FieldDescr * field_descr,
 int monitor_iter, 
 int rank,
 int iter_max, double res_tol,
 bool diag_precon) 
  : Solver(monitor_iter),
    A_(NULL),
    M_((diag_precon) ?
       (Matrix *)(new EnzoMatrixDiagonal) :
       (Matrix *)(new EnzoMatrixIdentity)),
    rank_(rank),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    rr0_(0),
    rr_min_(0),rr_max_(0),
    ix_(0),  ib_(0),
    ir_(0), id_(0), iy_(0), iz_(0),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    rr_(0.0), rz_(0.0), rz2_(0.0), dy_(0.0), bs_(0.0), rs_(0.0), xs_(0.0),
    bc_(0.0)
{
  TRACE_SOLVER("EnzoSolverCg() ENTER");

  id_ = field_descr->field_id("D");
  ir_ = field_descr->field_id("R");
  iy_ = field_descr->field_id("Y");
  iz_ = field_descr->field_id("Z");

  /// Initialize default Refresh

  const int num_fields = field_descr->field_count();

  field_descr->ghost_depth    (ib_,&gx_,&gy_,&gz_);

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields(num_fields);

  /// Initialize matvec Refresh

  TRACE_SOLVER("EnzoSolverCg() EXIT");
}

//----------------------------------------------------------------------

EnzoSolverCg::~EnzoSolverCg() throw ()
{
  if (A_) delete A_;
  A_ = NULL;
  if (M_) delete M_;
  M_ = NULL;
}

//----------------------------------------------------------------------

void EnzoSolverCg::pup (PUP::er &p)
{
  TRACEPUP;

  Solver::pup(p);

  p | A_;
  p | M_;
  p | rank_;
  p | iter_max_;
  p | res_tol_;
  p | rr0_;
  p | rr_min_;
  p | rr_max_;
  p | ib_;
  p | ix_;
  p | ir_;
  p | id_;
  p | iy_;
  p | iz_;

  p | nx_;
  p | ny_;
  p | nz_;

  p | mx_;
  p | my_;
  p | mz_;

  p | gx_;
  p | gy_;
  p | gz_;

  p | iter_;
  
  p | rz_;
  p | rz2_;
  p | dy_;
  p | bs_;
  p | bc_;
  p | id_refresh_matvec_;

}

//======================================================================

void EnzoSolverCg::apply ( Matrix * A, int ix, int ib, Block * block) throw()
{
  Solver::begin_(block);

  A_ = A;
  ix_ = ix;
  ib_ = ib;
  
  TRACE_SOLVER("compute ENTER");
  Field field = block->data()->field();

  field.size           (&nx_,&ny_,&nz_);
  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  // assumes all fields involved in calculation have same precision
  int precision = field.precision(ib_);

  if      (precision == precision_single)
    compute_<float>      (enzo_block);
  else if (precision == precision_double)
    compute_<double>     (enzo_block);
  else if (precision == precision_quadruple)
    compute_<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverCg()", "precision %d not recognized", precision);
  TRACE_SOLVER("compute EXIT");

}

//======================================================================

extern CkReduction::reducerType sum_long_double_type;
extern CkReduction::reducerType sum_long_double_2_type;
extern CkReduction::reducerType sum_long_double_3_type;
extern CkReduction::reducerType sum_long_double_4_type;

template <class T>
void EnzoSolverCg::compute_ (EnzoBlock * enzo_block) throw()
//     X = initial guess
//     B = right-hand side
//     R = B - A*X
//     solve(M*Z = R)
//     D = Z
//     shift (B)
{
  TRACE_SOLVER("compute_() ENTER");

  iter_ = 0;

  Data * data = enzo_block->data();
  Field field = data->field();

  if (enzo_block->is_leaf()) {

    ///   - X = 0
    ///   - R = P = B ( residual with X = 0);

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
	  R[i] = B[i];
	}
      }
    }

    M_->matvec(id_,ir_,enzo_block);
    M_->matvec(iz_,ir_,enzo_block);
  }

  long double reduce[3] = {0.0};

  if (enzo_block->is_leaf()) {

    T * B = (T*) field.values(ib_);
    T * R = (T*) field.values(ir_);

    const int i0 = gx_ + mx_*(gy_ + my_*gz_);

    // reduce[0] = field.dot(ir_,ir_);
    // reduce[1] = sum_(B);
    // reduce[2] = count_();
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + (ix + mx_*(iy + my_*iz));
	  reduce[0] += R[i]*R[i];
	  reduce[1] += B[i];
	  ++reduce[2];
	}
      }
    }
    
  }

#ifdef DEBUG_SOLVER
  printf ("%s:%d %s DEBUG_SOLVER calling callback\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif
  TRACE_SOLVER("   loop_0a callback");
  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_0a<T>(NULL), 
		      enzo_block->proxy_array());
  TRACE_SOLVER("   loop_0a callback");

#ifdef DEBUG_SOLVER
  printf ("%s:%d %s DEBUG_SOLVER calling contribute\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif
	  
  TRACE_SOLVER("   loop_0a contribute");

  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  sum_long_double_3_type, 
			  callback);
  TRACE_SOLVER("   loop_0a contribute");
  
  TRACE_SOLVER("compute_() EXIT");
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_solver_cg_loop_0a (CkReductionMsg * msg)
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> refresh P for AP = MATVEC (A,P)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  TRACE_SOLVER("solver_cg_loop_0a ENTER");

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  long double * data = (long double *) msg->getData();

  solver->set_rr( data[0] );
  solver->set_bs( data[1] );
  solver->set_bc( data[2] );

  delete msg;
  
  // Refresh field faces then call r_solver_cg_matvec

  Refresh refresh (4,0,neighbor_leaf, sync_barrier);
  refresh.set_active(is_leaf());
  refresh.add_all_fields(this->data()->field().field_count());
  refresh_enter(CkIndex_EnzoBlock::r_solver_cg_matvec(NULL),&refresh);
  
  TRACE_SOLVER("solver_cg_loop_0a EXIT");
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_solver_cg_loop_0b (CkReductionMsg * msg)
/// ==> refresh P for AP = MATVEC (A,P)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  TRACE_SOLVER("solver_cg_loop_0b ENTER");

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  solver->set_iter ( ((int*)msg->getData())[0] );

  // Refresh field faces then call solver_matvec

  Refresh refresh (4,0,neighbor_leaf, sync_barrier);
  refresh.set_active(is_leaf());
  refresh.add_all_fields(this->data()->field().field_count());
  refresh_enter(CkIndex_EnzoBlock::r_solver_cg_matvec(NULL), &refresh);

  TRACE_SOLVER("solver_cg_loop_0b EXIT");
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_cg_matvec(CkReductionMsg * msg)
{
  delete msg;
  
  TRACE_SOLVER("solver_cg_matvec ENTER");
  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  Field field = data()->field();

  // assumes all fields involved in calculation have same precision
  int precision = field.precision(0);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (this);

  if      (precision == precision_single)    
    solver->shift_1<float>(enzo_block);
  else if (precision == precision_double)    
    solver->shift_1<double>(enzo_block);
  else if (precision == precision_quadruple) 
    solver->shift_1<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverCg()", "precision %d not recognized", precision);
  TRACE_SOLVER("solver_cg_matvec EXIT");
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverCg::shift_1 (EnzoBlock * enzo_block) throw()
{

  TRACE_SOLVER("Solver shift_1() ENTER");
  Data * data = enzo_block->data();
  Field field = data->field();

  if (enzo_block->is_leaf()) {

    cello::check(rr_,"rr_",__FILE__,__LINE__);
    cello::check(bs_,"bs_",__FILE__,__LINE__);
    cello::check(bc_,"bc_",__FILE__,__LINE__);

    T * B  = (T*) field.values(ib_);
    T * R  = (T*) field.values(ir_);

    if (iter_ == 0 && A_->is_singular())  {

      // shift rhs B by projection of B onto e: B~ <== B - (e*eT)/(eT*e) b
      // eT*e == n === zone count (bc)
      // eT*b == sum_i=1,n B[i]

      // shift_ (R,shift,R);
      // shift_ (B,shift,B);
      T shift = -bs_ / bc_;
      for (int iz=0; iz<mz_; iz++) {
	for (int iy=0; iy<my_; iy++) {
	  for (int ix=0; ix<mx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    R[i] += shift;
	    B[i] += shift;
	  }
	}
      }

      M_->matvec(id_,ir_,enzo_block);
      M_->matvec(iz_,ir_,enzo_block);
    } 
  }

  long double reduce = 0;

  if (enzo_block->is_leaf()) {

    T * R  = (T*) field.values(ir_);
    // reduce = field.dot(ir_,ir_);
    const int i0 = gx_ + mx_*(gy_ + my_*gz_);
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + (ix + mx_*(iy + my_*iz));
	  reduce += R[i]*R[i];
	}
      }
    }
  } 

  TRACE_SOLVER("   shift_1 callback");
  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_shift_1<T>(NULL), 
		      enzo_block->proxy_array());
  TRACE_SOLVER("   shift_1 callback");

#ifdef DEBUG_SOLVER
  printf ("%s:%d %s DEBUG_SOLVER calling contribute\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif

  TRACE_SOLVER("   shift_1 contribute");

  enzo_block->contribute (sizeof(long double), &reduce, 
			  sum_long_double_type, 
			  callback);
  TRACE_SOLVER("   shift_1 contribute");
  TRACE_SOLVER("Solver::shift_1() EXIT");
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_solver_cg_shift_1 (CkReductionMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  TRACE_SOLVER("Block solver_cg_shift_1 ENTER");
  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  solver->set_rr( ((long double*)msg->getData())[0] );

  delete msg;

  solver -> loop_2a<T>(this);

  TRACE_SOLVER("Block solver_cg_shift_1 EXIT");
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverCg::loop_2a (EnzoBlock * enzo_block) throw()
{
  Refresh refresh (4,0,neighbor_leaf, sync_barrier);
  refresh.set_active(enzo_block->is_leaf());
  refresh.add_all_fields(enzo_block->data()->field().field_count());
  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_cg_loop_2<T>(NULL),&refresh);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_cg_loop_2 (CkReductionMsg * msg)
{
  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  solver->loop_2b<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverCg::loop_2b (EnzoBlock * enzo_block) throw()
{
  TRACE_SOLVER("Block loop_2b ENTER");
  cello::check(rr_,"rr_",__FILE__,__LINE__);

  if (iter_ == 0) {
    rr0_ = rr_;
    rr_min_ = rr_;
    rr_max_ = rr_;
  } else {
    rr_min_ = std::min(rr_min_,rr_);
    rr_max_ = std::max(rr_max_,rr_);
  }

  const bool is_root = enzo_block->index().is_root();

  if (is_root) {

    Monitor * monitor = enzo_block->simulation()->monitor();

    if (iter_ == 0) {
      monitor->print ("Enzo", "CG iter %04d  rr0 %g",
		      iter_,(double)(rr0_));
    }

    if (monitor_iter_ && (iter_ % monitor_iter_) == 0 ) {
      monitor_output_ (enzo_block,iter_,rr0_,rr_min_,rr_,rr_max_);
    }
  }

  if (rr_ / rr0_ < res_tol_) {

    if (is_root && monitor_iter_) {
      monitor_output_ (enzo_block,iter_,rr0_,rr_min_,rr_,rr_max_);
    }
    end<T> (enzo_block,return_converged);

  } else if (iter_ >= iter_max_)  {

    end<T>(enzo_block,return_error_max_iter_reached);

  } else {
    
    Data * data = enzo_block->data();
    Field field = data->field();

    if (enzo_block->is_leaf()) {

      double hx,hy,hz;
      data->field_cell_width(&hx,&hy,&hz);

      A_->matvec(iy_,id_,enzo_block);

    }

    long double reduce[3] = {0.0} ;

    if (enzo_block->is_leaf()) {

      T * D = (T*) field.values(id_);
      T * Y = (T*) field.values(iy_);
      T * R = (T*) field.values(ir_);
      T * Z = (T*) field.values(iz_);

      const int i0 = gx_ + mx_*(gy_ + my_*gz_);
	 
      for (int iz=0; iz<nz_; iz++) {
	for (int iy=0; iy<ny_; iy++) {
	  for (int ix=0; ix<nx_; ix++) {
	    int i = i0 + (ix + mx_*(iy + my_*iz));
	    reduce[0] += R[i]*R[i];
	    reduce[1] += R[i]*Z[i];
	    reduce[2] += D[i]*Y[i];
	  }
	}
      }
    }

    TRACE_SOLVER("loop_3 callback");

#ifdef DEBUG_SOLVER
    printf ("%s:%d %s DEBUG_SOLVER calling contribute\n",
	    __FILE__,__LINE__,enzo_block->name().c_str());
#endif

    TRACE_SOLVER("loop_3 contribute");
    CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_3<T>(NULL), 
			enzo_block->proxy_array());

    enzo_block->contribute (3*sizeof(long double), &reduce, 
			    sum_long_double_3_type,
			    callback);
    TRACE_SOLVER("loop_3 contribute");
  }
  TRACE_SOLVER("Block loop_2b EXIT");
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_solver_cg_loop_3 (CkReductionMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  TRACE_SOLVER("Block solver_cg_loop_3 ENTER");

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  long double * data = (long double *) msg->getData();

  solver->set_rr(data[0]);
  solver->set_rz(data[1]);
  solver->set_dy(data[2]);
  
  delete msg;

  solver -> loop_4<T>(this);

  TRACE_SOLVER("Block solver_cg_loop_3 EXIT");
  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverCg::loop_4 (EnzoBlock * enzo_block) throw ()
//  a = rz / dy;
//  X = X + a*D;
//  R = R - a*Y;
//  solve(M*Z = R)
//  rz2 = dot(R,Z)
//  b = rz2 / rz;
//  D = Z + b*D;
//  rz = rz2;
{
  TRACE_SOLVER("Solver loop_4 ENTER");

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

    //    field.axpy (ix_,  a ,id_, ix_);
    //    field.axpy (ir_, -a, iy_, ir_);
    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  X[i] += a * D[i];
	  R[i] -= a * Y[i];
	}
      }
    }

    M_->matvec(iz_,ir_,enzo_block);

  }

  long double reduce[3] = {0.0};

  if (enzo_block->is_leaf()) {

    T * X = (T*) field.values(ix_);
    T * R = (T*) field.values(ir_);
    T * Z = (T*) field.values(iz_);

    const int i0 = gx_ + mx_*(gy_ + my_*gz_);
       
    //    reduce[0] = field.dot(ir_,iz_);
    //    reduce[1] = sum_(R);
    //    reduce[2] = sum_(X);
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + (ix + mx_*(iy + my_*iz));
	  reduce[0] += R[i]*Z[i];
	  reduce[1] += R[i];
	  reduce[2] += X[i];
	}
      }
    }
  }

  TRACE_SOLVER("loop_5 callback");
  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_5<T>(NULL), 
		      enzo_block->proxy_array());
  TRACE_SOLVER("loop_5 callback");

#ifdef DEBUG_SOLVER
  printf ("%s:%d %s DEBUG_SOLVER calling contribute\n",
	  __FILE__,__LINE__,enzo_block->name().c_str());
#endif

  TRACE_SOLVER("loop_5 contribute");

  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  sum_long_double_3_type, 
			  callback);
  TRACE_SOLVER("loop_5 contribute");

  TRACE_SOLVER("Solver loop_4 EXIT");
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::r_solver_cg_loop_5 (CkReductionMsg * msg)
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> solver_cg_loop_6
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  TRACE_SOLVER("Block solver_cg_loop_5 ENTER");

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  long double * data = (long double *) msg->getData();

  solver->set_rz2(data[0]);
  solver->set_rs (data[1]);
  solver->set_xs (data[2]);

  delete msg;

  solver -> loop_6<T>(this);

  TRACE_SOLVER("Block solver_cg_loop_5 EXIT");

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverCg::loop_6 (EnzoBlock * enzo_block) throw ()
//  rz2 = dot(R,Z)
//  b = rz2 / rz;
//  D = Z + b*D;
//  rz = rz2;
{
  TRACE_SOLVER("Solver loop_6 ENTER");

  cello::check(rz2_,"rz2_",__FILE__,__LINE__);
  cello::check(rs_,"rs_",__FILE__,__LINE__);
  cello::check(xs_,"xs_",__FILE__,__LINE__);

  Field field = enzo_block->data()->field();

  if (enzo_block->is_leaf()) {

    if (A_->is_singular())  {

      // shift rhs B by projection of B onto e: B~ <== B - (e*eT)/(eT*e) b
      // eT*e == n === zone count (bc)
      // eT*b == sum_i=1,n B[i]

      T * X  = (T*) field.values(ix_);
      T * R  = (T*) field.values(ir_);

      // shift_ (X,T(-xs_/bc_),X);
      // shift_ (R,T(-rs_/bc_),R);

      for (int iz=0; iz<mz_; iz++) {
	for (int iy=0; iy<my_; iy++) {
	  for (int ix=0; ix<mx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    X[i] -= T(xs_/bc_);
	    R[i] -= T(rs_/bc_);
	  }
	}
      }
      
    }

    T * D  = (T*) field.values(id_);
    T * Z  = (T*) field.values(iz_);

    T b = rz2_ / rz_;

    cello::check(b,"b",__FILE__,__LINE__);

    // zaxpy_ (D,b,D,Z);
    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  D[i] = Z[i] + b * D[i];
	}
      }
    }
  }

  int iter = iter_ + 1;

  TRACE_SOLVER("loop_0b callback");
  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_0b<T>(NULL), 
		      enzo_block->proxy_array());
  TRACE_SOLVER("loop_0b callback");
    
  TRACE_SOLVER("loop_0b contribute");

  enzo_block->contribute (sizeof(int), &iter, 
			  CkReduction::max_int, callback);
  TRACE_SOLVER("loop_0b contribute");

  TRACE_SOLVER("Block loop_6 EXIT");
}

//----------------------------------------------------------------------


template <class T>
void EnzoSolverCg::end (EnzoBlock * enzo_block,int retval) throw ()
///    if (return == return_converged) {
///       ==> exit()
///    } else {
///       ERROR (return-)
///    }
{

  Solver::end_(enzo_block);
  
  TRACE_SOLVER("Solver end ENTER");

  CkCallback(callback_,
	     CkArrayIndexIndex(enzo_block->index()),
	     enzo_block->proxy_array()).send();

  TRACE_SOLVER("Solver end EXIT");
}

//----------------------------------------------------------------------

void EnzoSolverCg::exit_() throw()
/// deallocate temporary vectors
{
  TRACE_SOLVER("exit_");
}
