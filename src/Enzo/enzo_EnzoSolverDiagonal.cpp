// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDiagonal.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-03-01
/// @brief    Diagonal preconditioner: return X = B\diag(A)

#include "enzo.hpp"

#include "enzo.decl.h"

// #define DEBUG_SOLVER

#ifdef DEBUG_SOLVER
#   define TRACE_SOLVER(solver)						\
  CkPrintf ("%d %s:%d TRACE_DIAGONAL %s %p\n",CkMyPe(),__FILE__,__LINE__,solver,this); \
  fflush(stdout);
#else
#   define TRACE_SOLVER(solver) /*  */ 
#endif

EnzoSolverDiagonal::EnzoSolverDiagonal
(std::string name,
 std::string field_x,
 std::string field_b,
 int monitor_iter,
 int restart_cycle,
 int solve_type,
 int index_prolong,
 int index_restrict) throw()
  : Solver
    (name,
     field_x,
     field_b,
     monitor_iter,
     restart_cycle,
     solve_type,
     index_prolong,
     index_restrict)
{
  id_ = cello::field_descr()->insert_temporary();
}

//======================================================================

void EnzoSolverDiagonal::apply (std::shared_ptr<Matrix> A, Block * block) throw()
{
  Solver::begin_(block);

  TRACE_SOLVER("apply()");
  if (is_finest_(block)) {

    TRACE_SOLVER("is_finest");
    Field field = block->data()->field();

    // // assumes all fields involved in calculation have same precision
    // int precision = field.precision(ib);

    compute_(A,block);
  }
  
  Solver::end_(block);
}

//======================================================================

void EnzoSolverDiagonal::compute_
( std::shared_ptr<Matrix> A, Block * block) throw()
//     X = B / diag(A)
{
  TRACE_SOLVER("compute_() ENTER");

  Field field = block->data()->field();

  int mx,my,mz;
  field.dimensions (ib_,&mx,&my,&mz);

  if (is_finest_(block)) {

    field.allocate_temporary(id_);

    ///   - X = 0
    ///   - R = P = B ( residual with X = 0);

    A->diagonal(id_,block);

    enzo_float * X = (enzo_float*) field.values(ix_);
    enzo_float * B = (enzo_float*) field.values(ib_);
    enzo_float * D = (enzo_float*) field.values(id_);

    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  X[i] = B[i] / D[i];
	}
      }
    }
    field.deallocate_temporary(id_);
  }
  TRACE_SOLVER("compute_() EXIT");
}

