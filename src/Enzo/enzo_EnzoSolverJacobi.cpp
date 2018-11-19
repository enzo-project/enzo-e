// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverJacobi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoSolverJacobi class

#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_COPY
// #define DEBUG_SOLVER
// #define DEBUG_TRACE

#ifdef DEBUG_TRACE
#  define TRACE_JACOBI(BLOCK,METHOD) \
  CkPrintf ("%s:%d %s TRACE_JACOBI %s\n", \
	    __FILE__,__LINE__,BLOCK->name().c_str(),METHOD);
#else
#  define TRACE_JACOBI(BLOCK,METHOD) /* empty */
#endif


#ifdef DEBUG_SOLVER
#   define DEBUG_FIELD(BLOCK,IX,NAME)			\
  {									\
    Field field = BLOCK->data()->field();				\
    int mx,my,mz;							\
    field.dimensions(IX,&mx,&my,&mz);					\
    int gx,gy,gz;							\
    field.ghost_depth(IX,&gx,&gy,&gz);					\
    enzo_float * X = (enzo_float*) field.values(IX);			\
    double xx=0.0;							\
    double yy=0.0;							\
    for (int i=0; i<mx*my*mz; i++) {					\
      xx+=X[i]*X[i];							\
    }									\
    for (int iz=gz; iz<mz-gz; iz++) {					\
      for (int iy=gy; iy<my-gy; iy++) {					\
	for (int ix=gx; ix<mx-gx; ix++) {				\
	  int i = ix + mx*(iy + my*iz);					\
	    yy+=X[i]*X[i];						\
	  }								\
	}								\
      }									\
      CkPrintf ("%-8s %-10s %d DEBUG_FIELD ||%s|| = [%g] (%g)\n",		\
		BLOCK->name().c_str(),this->name().c_str(),	\
		__LINE__,NAME,xx,yy);				\
  }
#else
#   define DEBUG_FIELD(BLOCK,IX,NAME) /* ... */
#endif

//----------------------------------------------------------------------

EnzoSolverJacobi::EnzoSolverJacobi
( std::string name,
  std::string field_x,
  std::string field_b,
  int monitor_iter,
  int restart_cycle,
  int solve_type,
  double weight, int iter_max) throw()
  : Solver(name,
	   field_x,
	   field_b,
	   monitor_iter,
	   restart_cycle,
	   solve_type),
    A_ (NULL),
    ir_ (-1),
    id_ (-1),
    w_(weight),
    n_(iter_max)
{
  // Reserve temporary fields

  FieldDescr * field_descr = cello::field_descr();

  id_ = field_descr->insert_temporary();
  ir_ = field_descr->insert_temporary();

  ScalarDescr * scalar_descr_int = cello::scalar_descr_int();
  i_iter_ = scalar_descr_int->new_value(name_ + ":iter");  
}

//----------------------------------------------------------------------

void EnzoSolverJacobi::apply
( std::shared_ptr<Matrix> A, Block * block) throw()
{
  TRACE_JACOBI(block,"apply()");

  begin_(block);

  if (! is_finest_(block)) {
    // Exit if block does not participate in the solve
    TRACE_JACOBI(block,"end()");
    Solver::end_(block);

  } else {
 
    A_ = A;

    Field field = block->data()->field();

    allocate_temporary_(field,block);

    (*piter_(block)) = 0.0;
  
    // Refresh X
    Refresh refresh (4,0,neighbor_type_(), sync_type_(), sync_id_());

    refresh.add_field (ix_);
    
    block->refresh_enter
      (CkIndex_EnzoBlock::p_solver_jacobi_continue(),&refresh);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_jacobi_continue()
{
  TRACE_JACOBI(this,"p_solver_jacobi_continue()");
  
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverJacobi * solver = 
    static_cast<EnzoSolverJacobi *> (this->solver());

  solver->compute(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverJacobi::compute(Block * block)
{
  TRACE_JACOBI(block,"compute()");

  Field field = block->data()->field();

  if (*piter_(block) < n_) {

    apply_(block);

  } else {

    deallocate_temporary_ (field,block);

    Solver::end_(block);

  }
}  

//----------------------------------------------------------------------

void EnzoSolverJacobi::apply_(Block * block)
{
  TRACE_JACOBI(block,"apply_()");
  
  Field field = block->data()->field();

  int mx,my,mz;
  field.dimensions(ix_,&mx,&my,&mz);
  // int gx,gy,gz;
  // field.ghost_depth(ix_,&gx,&gy,&gz);

  const int ng = A_->ghost_depth();
  const int gx = (mx > 1) ? ng : 0;
  const int gy = (my > 1) ? ng : 0;
  const int gz = (mz > 1) ? ng : 0;

  A_->diagonal (id_, block,ng);
  A_->residual (ir_, ib_, ix_, block,ng);

#ifdef DEBUG_COPY
    {
      enzo_float * R = (enzo_float*) field.values(ir_);
      enzo_float * R_J = (enzo_float*) field.values("R_J");
      enzo_float * D = (enzo_float*) field.values(id_);
      enzo_float * D_J = (enzo_float*) field.values("D_J");
      enzo_float * X = (enzo_float*) field.values(ix_);
      enzo_float * X_J = (enzo_float*) field.values("X_J");
      enzo_float * B = (enzo_float*) field.values(ib_);
      enzo_float * B_J = (enzo_float*) field.values("B_J");
      double rsum=0.0;
      double dsum=0.0;
      double xsum=0.0;
      double bsum=0.0;
      for (int i=0; i<mx*my*mz; i++) {
	R_J[i]=R[i];
	D_J[i]=D[i];
	X_J[i]=X[i];
	B_J[i]=B[i];
	rsum+=std::abs(R[i]);
	dsum+=std::abs(D[i]);
	xsum+=X[i];
	bsum+=std::abs(B[i]);
      }
      CkPrintf ("DEBUG_COPY rsum dsum xsum bsum %g %g %g %g\n",rsum,dsum,xsum,bsum);
    }
#endif    
  enzo_float * X = (enzo_float*) field.values(ix_);
  enzo_float * R = (enzo_float*) field.values(ir_);
  enzo_float * D = (enzo_float*) field.values(id_);

  if (w_ == 1.0) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  X[i] += R[i] / D[i];
	}
      }
    }
  } else {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  X[i] = w_*(R[i] / D[i]) + (1.0-w_)*X[i];
	}
      }
    }
  }

  // Next iteration

  (*piter_(block))++;
  
  // Refresh X
  Refresh refresh (4,0,neighbor_type_(), sync_type_(), sync_id_());

  refresh.add_field (ix_);

  block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_jacobi_continue(),&refresh);


}

