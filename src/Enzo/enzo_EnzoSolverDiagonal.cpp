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
( Matrix * A, int ix, int ib, Block * block) throw()
{
  Solver::begin_(block);

  TRACE_SOLVER("apply()");
  if (is_active_(block)) {

    TRACE_SOLVER("is_active");
    Field field = block->data()->field();

    // assumes all fields involved in calculation have same precision
    int precision = field.precision(ib);

    if      (precision == precision_single)
      compute_<float>      (A,ix,ib,block);
    else if (precision == precision_double)
      compute_<double>     (A,ix,ib,block);
    else if (precision == precision_quadruple)
      compute_<long double>(A,ix,ib,block);
    else 
      ERROR1("EnzoSolverDiagonal()", "precision %d not recognized", precision);
    TRACE_SOLVER("compute EXIT");
  }
  
  Solver::end_(block);
  CkCallback(callback_,
	     CkArrayIndexIndex(block->index()),
	     block->proxy_array()).send();

}

//======================================================================

template <class T>
void EnzoSolverDiagonal::compute_
( Matrix * A, int ix, int ib, Block * block) throw()
//     X = B / diag(A)
{
  TRACE_SOLVER("compute_() ENTER");

  Field field = block->data()->field();

  int mx,my,mz;
  field.dimensions (ib,&mx,&my,&mz);


  if (is_active_(block)) {

    ///   - X = 0
    ///   - R = P = B ( residual with X = 0);

    const int id = field.field_id("diagonal_D");

    A->diagonal(id,block);

    T * X = (T*) field.values(ix);
    T * B = (T*) field.values(ib);
    T * D = (T*) field.values(id);

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

