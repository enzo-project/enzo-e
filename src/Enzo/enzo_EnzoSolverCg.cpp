// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-1-08
/// @brief    Implements the CG Krylov iterative linear solver

#include "enzo.hpp"

#include "enzo.decl.h"
// #define DEBUG_RESID
//----------------------------------------------------------------------

EnzoSolverCg::EnzoSolverCg 
(FieldDescr * field_descr,
 int monitor_iter, 
 int rank,
 int iter_max, double res_tol,
 int min_level, int max_level,
 int index_precon,
 bool local
 )
  : Solver(monitor_iter,min_level,max_level), 
    A_(NULL),
    ix_(0),  ib_(0),
    index_precon_(index_precon),
    rank_(rank),
    iter_max_(iter_max), 
    ir_(0), id_(0), iy_(0), iz_(0),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    res_tol_(res_tol),
    rr0_(0.0),
    rr_min_(0.0),rr_max_(0.0),
    rr_(0.0), rz_(0.0), rz2_(0.0), dy_(0.0), bs_(0.0), rs_(0.0), xs_(0.0),
    bc_(0.0),
    local_(local)
{
  id_ = field_descr->insert_temporary();
  ir_ = field_descr->insert_temporary();
  iy_ = field_descr->insert_temporary();
  iz_ = field_descr->insert_temporary();

  /// Initialize default Refresh

  field_descr->ghost_depth    (ib_,&gx_,&gy_,&gz_);

  if (! local_) {
    const int ir = add_refresh(4,0,neighbor_type_(),
			       sync_type_(),
			       enzo_sync_id_solver_cg);
    if (CkMyPe() == 0) {
      WARNING("EnzoSolverCg",
	      "assuming solution field is 'potential'");
    }
    refresh(ir)->add_field (field_descr->field_id("potential"));
    
    refresh(ir)->add_field (id_);
    refresh(ir)->add_field (ir_);
    refresh(ir)->add_field (iy_);
    refresh(ir)->add_field (iz_);
  }

}

//----------------------------------------------------------------------

EnzoSolverCg::~EnzoSolverCg() throw ()
{
  if (A_) delete A_;
  A_ = NULL;
}

//----------------------------------------------------------------------

void EnzoSolverCg::pup (PUP::er &p)
{
  TRACEPUP;

  Solver::pup(p);

  p | A_;
  p | ix_;
  p | ib_;

  p | index_precon_;
  
  p | rank_;
  p | iter_max_;
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
  
  p | res_tol_;
  p | rr0_;
  p | rr_min_;
  p | rr_max_;
  p | rz_;
  p | rz2_;
  p | dy_;
  p | bs_;
  p | bc_;

  p | local_;
}

//======================================================================

void EnzoSolverCg::apply ( Matrix * A, int ix, int ib, Block * block) throw()
{
  Solver::begin_(block);

  A_ = A;
  ix_ = ix;
  ib_ = ib;
  
  Field field = block->data()->field();

  allocate_temporary_(field,block);

  field.size           (&nx_,&ny_,&nz_);
  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  // assumes all fields involved in calculation have same precision
  // int precision = field.precision(ib_);

  compute_(enzo_block);

}

//======================================================================

void EnzoSolverCg::compute_ (EnzoBlock * enzo_block) throw()
//     X = initial guess
//     B = right-hand side
//     R = B - A*X
//     solve(M*Z = R)
//     D = Z
//     shift (B)
{

  // If local, call serial CG solver
  if (local_) {
    local_cg_(enzo_block);
    return;
  }
  
  iter_ = 0;

  Field field = enzo_block->data()->field();

  if (is_active_(enzo_block)) {

    ///   - X = 0
    ///   - R = P = B ( residual with X = 0);

    enzo_float * B = (enzo_float*) field.values(ib_);
    enzo_float * X = (enzo_float*) field.values(ix_);
    enzo_float * R = (enzo_float*) field.values(ir_);
    enzo_float * D = (enzo_float*) field.values(id_);
    enzo_float * Z = (enzo_float*) field.values(iz_);

    for (int i=0; i<mx_*my_*mz_; i++) {
      X[i] = 0.0;
      R[i] = B[i];
      D[i] = R[i];
      Z[i] = R[i];
    }

  }

  long double reduce[3] = {0.0};

  if (is_active_(enzo_block)) {

    enzo_float * B = (enzo_float*) field.values(ib_);
    enzo_float * R = (enzo_float*) field.values(ir_);

    const int i0 = gx_ + mx_*(gy_ + my_*gz_);

    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + (ix + mx_*(iy + my_*iz));
	  reduce[0] += R[i]*R[i];
	  reduce[1] += B[i];
	}
      }
    }
    reduce[2] = nx_*ny_*nz_;
  }

  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_0a(NULL), 
		      enzo_block->proxy_array());
	  
  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  sum_long_double_3_type, 
			  callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_cg_loop_0a (CkReductionMsg * msg)
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> refresh P for AP = MATVEC (A,P)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  solver->loop_0a(this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverCg::loop_0a
(EnzoBlock * enzo_block, CkReductionMsg * msg) throw ()
{

  long double * data = (long double *) msg->getData();

  rr_ = data[0];
  bs_ = data[1];
  bc_ = data[2];

  delete msg;

// Refresh field faces then call r_solver_cg_matvec

  Refresh refresh (4,0,neighbor_type_(), sync_type_(),
		   enzo_sync_id_solver_cg_loop_0a);
  refresh.set_active(is_active_(enzo_block));

  //refresh.add_all_fields();

  refresh.add_field (id_);
  refresh.add_field (ir_);
  refresh.add_field (iy_);
  refresh.add_field (iz_);
  
  enzo_block->refresh_enter(CkIndex_EnzoBlock::r_solver_cg_matvec(),&refresh);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_cg_loop_0b (CkReductionMsg * msg)
/// ==> refresh P for AP = MATVEC (A,P)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  solver->loop_0b(this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverCg::loop_0b
(EnzoBlock * enzo_block, CkReductionMsg * msg) throw ()
{

  set_iter ( ((int*)msg->getData())[0] );

  delete msg;
  
  // Refresh field faces then call solver_matvec

  Refresh refresh (4,0,neighbor_type_(), sync_type_(),
		   enzo_sync_id_solver_cg_loop_0b);
  refresh.set_active(is_active_(enzo_block));
  refresh.add_all_fields();
  refresh.add_field (id_);
  refresh.add_field (ir_);
  refresh.add_field (iy_);
  refresh.add_field (iz_);

  enzo_block->refresh_enter(CkIndex_EnzoBlock::r_solver_cg_matvec(), &refresh);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_cg_matvec()
{
  
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  //  Field field = data()->field();

  // assumes all fields involved in calculation have same precision
  //  int precision = field.precision(0);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (this);

  solver->shift_1(enzo_block);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverCg::shift_1 (EnzoBlock * enzo_block) throw()
{
  Data * data = enzo_block->data();
  Field field = data->field();

  if (is_active_(enzo_block)) {

    cello::check(rr_,"rr_",__FILE__,__LINE__);
    cello::check(bs_,"bs_",__FILE__,__LINE__);
    cello::check(bc_,"bc_",__FILE__,__LINE__);

    enzo_float * B  = (enzo_float*) field.values(ib_);
    enzo_float * R  = (enzo_float*) field.values(ir_);

    if (iter_ == 0 && A_->is_singular())  {

      // shift rhs B by projection of B onto e: B~ <== B - (e*eT)/(eT*e) b
      // eT*e == n === zone count (bc)
      // eT*b == sum_i=1,n B[i]

      // shift_ (R,shift,R);
      // shift_ (B,shift,B);
      cello::check(bc_,"bc_",__FILE__,__LINE__);
      cello::check(bs_,"bs_",__FILE__,__LINE__);
  
      enzo_float shift = -bs_ / bc_;
      cello::check(shift,"shift",__FILE__,__LINE__);
      enzo_float * D = (enzo_float*) field.values(id_);
      enzo_float * Z = (enzo_float*) field.values(iz_);
      for (int i=0; i<mx_*my_*mz_; i++) {
	R[i] += shift;
	B[i] += shift;
	D[i] = R[i];
	Z[i] = R[i];
      }
    } 
  }

  long double reduce = 0;

  if (is_active_(enzo_block)) {

    enzo_float * R  = (enzo_float*) field.values(ir_);
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

  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_shift_1(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute (sizeof(long double), &reduce, 
			  sum_long_double_type, 
			  callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_cg_shift_1 (CkReductionMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  solver->set_rr( ((long double*)msg->getData())[0] );

  delete msg;

  solver -> loop_2a(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverCg::loop_2a (EnzoBlock * enzo_block) throw()
{
  Refresh refresh (4,0,neighbor_type_(), sync_type_(),
		   enzo_sync_id_solver_cg_loop_2a);
  refresh.set_active(is_active_(enzo_block));
  refresh.add_all_fields();
  refresh.add_field (id_);
  refresh.add_field (ir_);
  refresh.add_field (iy_);
  refresh.add_field (iz_);

  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_cg_loop_2(),&refresh);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_cg_loop_2 ()
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  
  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  solver->loop_2b(this);
  performance_stop_(perf_compute,__FILE__,__LINE__);
  
}

//----------------------------------------------------------------------

void EnzoSolverCg::loop_2b (EnzoBlock * enzo_block) throw()
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

  const bool is_root = enzo_block->index().is_root();
  if (is_root) monitor_output_(enzo_block);

  const bool is_converged = (rr_ / rr0_ < res_tol_);
  const bool is_diverged = (iter_ >= iter_max_);
    
  if (is_converged) {

    end (enzo_block,return_converged);

  } else if (is_diverged)  {

    end (enzo_block,return_error);

  } else {

    // else continue
    Data * data = enzo_block->data();
    Field field = data->field();

    if (is_active_(enzo_block)) {

      A_->matvec(iy_,id_,enzo_block);

    }

    long double reduce[3] = {0.0} ;

    if (is_active_(enzo_block)) {

      enzo_float * D = (enzo_float*) field.values(id_);
      enzo_float * Y = (enzo_float*) field.values(iy_);
      enzo_float * R = (enzo_float*) field.values(ir_);
      enzo_float * Z = (enzo_float*) field.values(iz_);

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

    CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_3(NULL), 
			enzo_block->proxy_array());

    enzo_block->contribute (3*sizeof(long double), &reduce, 
			    sum_long_double_3_type,
			    callback);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_cg_loop_3 (CkReductionMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  long double * data = (long double *) msg->getData();

  solver->set_rr(data[0]);
  solver->set_rz(data[1]);
  solver->set_dy(data[2]);
  
  delete msg;

  solver -> loop_4(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

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

  cello::check(rr_,"rr_",__FILE__,__LINE__);
  cello::check(rz_,"rz_",__FILE__,__LINE__);
  cello::check(dy_,"dy_",__FILE__,__LINE__);

  Data * data = enzo_block->data();
  Field field = data->field();

  if (is_active_(enzo_block)) {

    enzo_float * X = (enzo_float*) field.values(ix_);
    enzo_float * D = (enzo_float*) field.values(id_);
    enzo_float * R = (enzo_float*) field.values(ir_);
    enzo_float * Y = (enzo_float*) field.values(iy_);

    enzo_float a = rz_ / dy_;

    cello::check(a,"a",__FILE__,__LINE__);

    //    field.axpy (ix_,  a ,id_, ix_);
    //    field.axpy (ir_, -a, iy_, ir_);
    for (int i=0; i<mx_*my_*mz_; i++) {
      X[i] += a * D[i];
      R[i] -= a * Y[i];
    }

    enzo_float * Z = (enzo_float*) field.values(iz_);
    
    // M_->matvec(iz_,ir_,enzo_block);
    for (int i=0; i<mx_*my_*mz_; i++) {
      Z[i] = R[i];
    }

#ifdef DEBUG_RESID
    CkPrintf ("Copying residual %s\n",enzo_block->name().c_str());
    enzo_float * residual = (enzo_float*) field.values("residual");
    for (int iz=0; iz<nz_; iz++) {
      int kz=iz+gz_;
      for (int iy=0; iy<ny_; iy++) {
	int ky=iy+gy_;
	for (int ix=0; ix<nx_; ix++) {
	  int kx=ix+gx_;
	  int i = kx + mx_*(ky + my_*kz);
	  residual[i]=R[i];
	}
      }
    }
#endif    
  }

  long double reduce[3] = {0.0};

  if (is_active_(enzo_block)) {

    enzo_float * X = (enzo_float*) field.values(ix_);
    enzo_float * R = (enzo_float*) field.values(ir_);
    enzo_float * Z = (enzo_float*) field.values(iz_);

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

  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_5(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute (3*sizeof(long double), &reduce, 
			  sum_long_double_3_type, 
			  callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_cg_loop_5 (CkReductionMsg * msg)
/// - EnzoBlock accumulate global contribution to DOT(R,R)
/// ==> solver_cg_loop_6
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverCg * solver = 
    static_cast<EnzoSolverCg*> (this->solver());

  long double * data = (long double *) msg->getData();

  solver->set_rz2(data[0]);
  solver->set_rs (data[1]);
  solver->set_xs (data[2]);

  delete msg;

  solver -> loop_6(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverCg::loop_6 (EnzoBlock * enzo_block) throw ()
//  rz2 = dot(R,Z)
//  b = rz2 / rz;
//  D = Z + b*D;
//  rz = rz2;
{
  cello::check(rz2_,"rz2_",__FILE__,__LINE__);
  cello::check(rs_,"rs_",__FILE__,__LINE__);
  cello::check(xs_,"xs_",__FILE__,__LINE__);

  Field field = enzo_block->data()->field();

  if (is_active_(enzo_block)) {

    if (A_->is_singular())  {

      // shift rhs B by projection of B onto e: B~ <== B - (e*eT)/(eT*e) b
      // eT*e == n === zone count (bc)
      // eT*b == sum_i=1,n B[i]

      enzo_float * X  = (enzo_float*) field.values(ix_);
      enzo_float * R  = (enzo_float*) field.values(ir_);

      // shift_ (X,T(-xs_/bc_),X);
      // shift_ (R,T(-rs_/bc_),R);

      for (int i=0; i<mx_*my_*mz_; i++) {
	X[i] -= enzo_float(xs_/bc_);
	R[i] -= enzo_float(rs_/bc_);
      }
      
   }

    enzo_float * D  = (enzo_float*) field.values(id_);
    enzo_float * Z  = (enzo_float*) field.values(iz_);

    enzo_float b = rz2_ / rz_;

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

  CkCallback callback(CkIndex_EnzoBlock::r_solver_cg_loop_0b(NULL), 
		      enzo_block->proxy_array());

  enzo_block->contribute (sizeof(int), &iter, 
			  CkReduction::max_int, callback);
}

//----------------------------------------------------------------------

void EnzoSolverCg::local_cg_(EnzoBlock * enzo_block)
{
  if (!is_active_(enzo_block)) {
    end(enzo_block,return_unknown);
    return;
  }
  
  iter_ = 0;

  Field field = enzo_block->data()->field();

  enzo_float * B = (enzo_float*) field.values(ib_);
  enzo_float * D = (enzo_float*) field.values(id_);
  enzo_float * R = (enzo_float*) field.values(ir_);
  enzo_float * X = (enzo_float*) field.values(ix_);
  enzo_float * Y = (enzo_float*) field.values(iy_);
  enzo_float * Z = (enzo_float*) field.values(iz_);

  for (int i=0; i<mx_*my_*mz_; i++) {
    X[i] = 0.0;
    R[i] = B[i];
    D[i] = R[i];
    Z[i] = R[i];
  }
  bs_ = 0.0;
  bc_ = 0.0;

  refresh_local_(ib_,enzo_block);
  refresh_local_(ix_,enzo_block);
  refresh_local_(ir_,enzo_block);
  refresh_local_(id_,enzo_block);
  refresh_local_(iz_,enzo_block);

  const int i0 = gx_ + mx_*(gy_ + my_*gz_);

  // Compute shift and update B if needed
  if (iter_ == 0 && A_->is_singular()) {
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + (ix + mx_*(iy + my_*iz));
	  bs_ += B[i];
	}
      }
    }
    bc_ = nx_*ny_*nz_;
    cello::check(bc_,"bc_",__FILE__,__LINE__);
    cello::check(bs_,"bs_",__FILE__,__LINE__);
    enzo_float shift = -bs_ / bc_;
    cello::check(shift,"shift",__FILE__,__LINE__);
    for (int i=0; i<mx_*my_*mz_; i++) {
      R[i] += shift;
      B[i] += shift;
      D[i] = R[i];
      Z[i] = R[i];
    }
  }

  // compute residual
  rr_ = 0.0;
  for (int iz=0; iz<nz_; iz++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + (ix + mx_*(iy + my_*iz));
	rr_ += R[i]*R[i];
      }
    }
  }

  cello::check(rr_,"rr_",__FILE__,__LINE__);
  cello::check(bs_,"bs_",__FILE__,__LINE__);
  cello::check(bc_,"bc_",__FILE__,__LINE__);

  rr0_ = rr_;
  
  bool is_converged = (rr_ / rr0_ < res_tol_);
  bool is_diverged = iter_ >= iter_max_;

  while ( (! is_converged) && (! is_diverged) ) {

    rr_min_ = std::min(rr_min_,rr_);
    rr_max_ = std::max(rr_max_,rr_);

    refresh_local_(id_,enzo_block);

    A_->matvec(iy_,id_,enzo_block);

    rr_ = 0.0;
    rz_ = 0.0;
    dy_ = 0.0;
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + (ix + mx_*(iy + my_*iz));
	  rr_ += R[i]*R[i];
	  rz_ += R[i]*Z[i];
	  dy_ += D[i]*Y[i];
	}
      }
    }

    cello::check(rr_,"rr_",__FILE__,__LINE__);
    cello::check(rz_,"rz_",__FILE__,__LINE__);
    cello::check(dy_,"dy_",__FILE__,__LINE__);

    enzo_float a = rz_ / dy_;

    cello::check(a,"a",__FILE__,__LINE__);

    for (int i=0; i<mx_*my_*mz_; i++) {
      X[i] += a * D[i];
      R[i] -= a * Y[i];
      Z[i] = R[i];
    }

    rz2_ = 0.0;
    rs_ = 0.0;
    xs_ = 0.0;
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + (ix + mx_*(iy + my_*iz));
	  rz2_ += R[i]*Z[i];
	  rs_  += R[i];
	  xs_  += X[i];
	}
      }
    }

    cello::check(rz2_,"rz2_",__FILE__,__LINE__);
    cello::check(rs_,"rs_",__FILE__,__LINE__);
    cello::check(xs_,"xs_",__FILE__,__LINE__);
    
    if (A_->is_singular()) {
      for (int i=0; i<mx_*my_*mz_; i++) {
	X[i] -= enzo_float(xs_/bc_);
	R[i] -= enzo_float(rs_/bc_);
      }
    }


    enzo_float b = rz2_ / rz_;

    cello::check(b,"b",__FILE__,__LINE__);

    for (int i=0; i<mx_*my_*mz_; i++) {
      D[i] = Z[i] + b * D[i];
    }
    
    monitor_output_(enzo_block);

    ++iter_;
    is_converged = (rr_ / rr0_ < res_tol_);
    is_diverged = iter_ >= iter_max_;
  }

  if (is_converged) {

    end (enzo_block,return_converged);

  } else if (is_diverged)  {

    end(enzo_block,return_error);

  }
    
}

//----------------------------------------------------------------------

void EnzoSolverCg::refresh_local_(int ix,EnzoBlock * enzo_block)
{

  enzo_float * X = (enzo_float*) enzo_block->data()->field().values(ix);

  // ASSUMES SINGULAR MATRIX IMPLIES PERIODIC DOMAIN.

  if (A_->is_singular()) {

    // XM ghost <- XP face (y)(z)
    for (int iz=gz_; iz<nz_+gz_; iz++) {
      for (int iy=gy_; iy<ny_+gy_; iy++) {
	for (int ix=0; ix<gx_; ix++) {
	  int is = (ix+nx_) + mx_*(iy + my_*iz);
	  int id = (ix    ) + mx_*(iy + my_*iz);
	  X[id] = X[is];
	}
      }
    }

    // XM ghost <- XP face (y)(z)
    for (int iz=gz_; iz<nz_+gz_; iz++) {
      for (int iy=gy_; iy<ny_+gy_; iy++) {
	for (int ix=nx_+gx_; ix<nx_+2*gx_; ix++) {
	  int is = (ix-nx_) + mx_*(iy + my_*iz);
	  int id = (ix    ) + mx_*(iy + my_*iz);
	  X[id] = X[is];
	}
      }
    }

    // YM ghost <- YP face [x](z)
    for (int iz=gz_; iz<nz_+gz_; iz++) {
      for (int iy=0; iy<gy_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int is = ix + mx_*((iy+ny_) + my_*iz);
	  int id = ix + mx_*((iy    ) + my_*iz);
	  X[id] = X[is];
	}
      }
    }

    // YP ghost <- YM face [x](z)
    for (int iz=gz_; iz<nz_+gz_; iz++) {
      for (int iy=ny_+gy_; iy<ny_+2*gy_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int is = ix + mx_*((iy-ny_) + my_*iz);
	  int id = ix + mx_*((iy    ) + my_*iz);
	  X[id] = X[is];
	}
      }
    }

    // ZM ghost <- ZP face [x][y]
    for (int iz=0; iz<gz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int is = ix + mx_*(iy + my_*(iz+nz_));
	  int id = ix + mx_*(iy + my_*(iz    ));
	  X[id] = X[is];
	}
      }
    }

    // ZP ghost <- ZM face [x][y]
    for (int iz=nz_+gz_; iz<nz_+2*gz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int is = ix + mx_*(iy + my_*(iz-nz_));
	  int id = ix + mx_*(iy + my_*(iz    ));
	  X[id] = X[is];
	}
      }
    }


  } else {

    ERROR("EnzoSolverCg::refresh_local_()",
	  "Only periodic boundary conditions available");

  }
  
}

//----------------------------------------------------------------------

void EnzoSolverCg::end (EnzoBlock * enzo_block,int retval) throw ()
///    if (return == return_converged) {
///       ==> exit()
///    } else {
///       ERROR (return-)
///    }
{
  Field field = enzo_block->data()->field();

  deallocate_temporary_(field,enzo_block);

  Solver::end_(enzo_block);
  
  CkCallback(callback_,
	     CkArrayIndexIndex(enzo_block->index()),
	     enzo_block->proxy_array()).send();
}

//----------------------------------------------------------------------

void EnzoSolverCg::monitor_output_(EnzoBlock * enzo_block)
{
  const bool l_first_iter = (iter_ == 0);
  const bool l_max_iter   = (iter_ >= iter_max_);
  const bool l_monitor    = (monitor_iter_ && (iter_ % monitor_iter_) == 0 );
  const bool l_converged  = (rr_ / rr0_ < res_tol_);

  const bool l_output =
    l_first_iter || l_max_iter || l_monitor || l_converged;
      
  if (l_output) {
    Solver::monitor_output_ (enzo_block,iter_,rr0_,rr_min_,rr_,rr_max_);
  }

}
