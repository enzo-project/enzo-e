// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverJacobi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoSolverJacobi class

#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_COPY
// #define DEBUG_TRACE

#ifdef DEBUG_TRACE
#  define TRACE_JACOBI(BLOCK,METHOD) \
  CkPrintf ("%s:%d %s TRACE_JACOBI %s\n", \
	    __FILE__,__LINE__,BLOCK->name().c_str(),METHOD);
#else
#  define TRACE_JACOBI(BLOCK,METHOD) /* empty */
#endif

//----------------------------------------------------------------------

EnzoSolverJacobi::EnzoSolverJacobi
( FieldDescr * field_descr, double weight, int iter_max) throw()
  : Solver(0),
    A_ (NULL),
    ix_ (-1),
    ib_ (-1),
    ir_ (-1),
    id_ (-1),
    w_(weight),
    n_(iter_max)
{
  // Reserve temporary fields
  id_ = field_descr->insert_temporary();
  ir_ = field_descr->insert_temporary();
}

//----------------------------------------------------------------------

void EnzoSolverJacobi::apply
( Matrix * A, int ix, int ib, Block * block) throw()
{
  TRACE_JACOBI(block,"apply()");
  
  begin_(block);

  A_ = A;
  ix_ = ix;
  ib_ = ib;

  Field field = block->data()->field();

  allocate_temporary_(field,block);
  
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->jacobi_iter_clear();
    
  // Refresh X
  Refresh refresh (4,0,neighbor_type_(), sync_type_(), sync_id_());

  refresh.add_field (ix_);

  block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_jacobi_continue(),&refresh);
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
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = block->data()->field();

  if (enzo_block->jacobi_iter() < n_) {

    apply_(block);

  } else {

    deallocate_temporary_ (field,block);

    Solver::end_(block);

    CkCallback(callback_,
	       CkArrayIndexIndex(block->index()),
	       block->proxy_array()).send();

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

  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
	int i = ix + mx*(iy + my*iz);
	X[i] += R[i] / D[i];
      }
    }
  }

  // Next iteration

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->jacobi_iter_increment();
  
  // Refresh X
  Refresh refresh (4,0,neighbor_type_(), sync_type_(), sync_id_());

  refresh.add_field (ix_);

  block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_jacobi_continue(),&refresh);


}

