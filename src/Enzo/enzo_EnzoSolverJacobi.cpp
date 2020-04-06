// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverJacobi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoSolverJacobi class

#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_COPY
// #define DEBUG_SOLVER
// #define DEBUG_NEW_REFRESH
// #define DEBUG_TRACE
// #define DEBUG_TRACE_CYCLE 0

#ifdef DEBUG_TRACE
#  define TRACE_JACOBI(BLOCK,SOLVER,METHOD)			\
  if (BLOCK->cycle() >= DEBUG_TRACE_CYCLE) {			\
    CkPrintf ("%s:%d %s %s TRACE_JACOBI active %d %s\n",			\
	      __FILE__,__LINE__,BLOCK->name().c_str(),(SOLVER?SOLVER->name().c_str():"Unknown"), \
              (SOLVER?SOLVER->is_finest(BLOCK):-1),METHOD);             \
  }
#else
#  define TRACE_JACOBI(BLOCK,SOLVER,METHOD) /* empty */
#endif


#ifdef DEBUG_SOLVER
#   define DEBUG_FIELD(BLOCK,IX,NAME)					\
  {									\
    Field field = BLOCK->data()->field();		\
  int mx,my,mz;								\
  field.dimensions(IX,&mx,&my,&mz);					\
  int gx,gy,gz;								\
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
    n_(iter_max),
    ir_smooth_(-1)
{
  // Reserve temporary fields

  id_ = cello::field_descr()->insert_temporary();
  ir_ = cello::field_descr()->insert_temporary();


  Refresh * refresh = cello::refresh(ir_post_);
  cello::simulation()->new_refresh_set_name(ir_post_,name);
#ifdef DEBUG_NEW_REFRESH  
  CkPrintf ("DEBUG_NEW_REFRESH %s:%d ir_post=%d\n",__FILE__,__LINE__,ir_post_);
#endif  

  refresh->add_field (ix_);

  ScalarDescr * scalar_descr_int = cello::scalar_descr_int();
  i_iter_ = scalar_descr_int->new_value(name_ + ":iter");  


  ir_smooth_ = add_new_refresh_();

#ifdef DEBUG_NEW_REFRESH  
  CkPrintf ("DEBUG_NEW_REFRESH %s:%d ir_smooth=%d\n",__FILE__,__LINE__,ir_smooth_);
#endif  

  Refresh * refresh_smooth = cello::refresh(ir_smooth_);
  cello::simulation()->new_refresh_set_name(ir_smooth_,name+":smooth");
  
  refresh_smooth->add_field (ix_);
  refresh_smooth->set_solver_id(index());
#ifdef DEBUG_NEW_REFRESH  
  CkPrintf ("DEBUG_NEW_REFRESH %s:%d id_solver=%d\n",__FILE__,__LINE__,index());
#endif  
  refresh_smooth->set_callback(CkIndex_EnzoBlock::p_solver_jacobi_continue());

}

//----------------------------------------------------------------------

void EnzoSolverJacobi::apply
( std::shared_ptr<Matrix> A, Block * block) throw()
{
  TRACE_JACOBI(block,this,"apply()");

  begin_(block);

  if (solve_type_ == solve_level && ! is_finest_(block))
    Solver::end_(block);


  A_ = A;

  Field field = block->data()->field();

  allocate_temporary_(field,block);

  (*piter_(block)) = 0.0;
  
  // Refresh X

  do_refresh_(block);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_jacobi_continue()
{
 
  performance_start_(perf_compute,__FILE__,__LINE__);

  EnzoSolverJacobi * solver = nullptr;  
  TRACE_JACOBI(this,solver,"p_solver_jacobi_continue()");

  solver = static_cast<EnzoSolverJacobi *> (this->solver());
  TRACE_JACOBI(this,solver,"p_solver_jacobi_continue()");

  solver->compute(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverJacobi::compute(Block * block)
{
  TRACE_JACOBI(block,this,"compute()");

  if (*piter_(block) < n_) {

    apply_(block);

  } else {

    Field field = block->data()->field();
    deallocate_temporary_ (field,block);

    TRACE_JACOBI(block,this,"end()");
    Solver::end_(block);

  }
}  

//----------------------------------------------------------------------

void EnzoSolverJacobi::apply_(Block * block)
{
  TRACE_JACOBI(block,this,"apply_()");
  
  Field field = block->data()->field();

  int mx,my,mz;
  field.dimensions(ix_,&mx,&my,&mz);
  // int gx,gy,gz;
  // field.ghost_depth(ix_,&gx,&gy,&gz);

  const int ng = A_->ghost_depth();
  const int gx = (mx > 1) ? ng : 0;
  const int gy = (my > 1) ? ng : 0;
  const int gz = (mz > 1) ? ng : 0;

  if (is_finest_(block)) {

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
  }
  // Next iteration

  (*piter_(block))++;
  
  // Refresh X

  do_refresh_(block);

}

//----------------------------------------------------------------------

void EnzoSolverJacobi::do_refresh_(Block * block)
{
  TRACE_JACOBI(block,this,"do_refresh()");
#ifdef DEBUG_NEW_REFRESH  
  CkPrintf ("DEBUG_NEW_REFRESH %s:%d ir_smooth=%d\n",__FILE__,__LINE__,ir_smooth_);
#endif  
  Refresh * refresh = cello::refresh(ir_smooth_);

  refresh->set_active(is_finest_(block));
  refresh->add_field (ix_);
  
#ifdef DEBUG_NEW_REFRESH  
  CkPrintf ("DEBUG_NEW_REFRESH %s:%d ir_smooth=%d\n",__FILE__,__LINE__,ir_smooth_);
#endif  
  block->new_refresh_start (ir_smooth_,
			    CkIndex_EnzoBlock::p_solver_jacobi_continue());
}

//----------------------------------------------------------------------
