// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDiagonal.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-03-01
/// @brief    Diagonal preconditioner: return X = B\diag(A)

#include "enzo.hpp"

#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

// #define DEBUG_SOLVER

#ifdef DEBUG_SOLVER
#   define TRACE_SOLVER(solver)						\
  CkPrintf ("%d %s:%d TRACE_DIAGONAL %s %p\n",CkMyPe(),__FILE__,__LINE__,solver,this); \
  fflush(stdout);
#else
#   define TRACE_SOLVER(solver) /*  */ 
#endif

//======================================================================

void EnzoSolverDiagonal::apply
( std::shared_ptr<Matrix> A, Block * block) throw()
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

    ///   - X = 0
    ///   - R = P = B ( residual with X = 0);

    const int id = field.field_id("diagonal_D");

    A->diagonal(id,block);

    enzo_float * X = (enzo_float*) field.values(ix_);
    enzo_float * B = (enzo_float*) field.values(ib_);
    enzo_float * D = (enzo_float*) field.values(id);

    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  X[i] = B[i] / D[i];
	}
      }
    }
  }
  TRACE_SOLVER("compute_() EXIT");
}

