// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-1-08
/// @brief    Implements the CG Krylov iterative linear solver

#include "enzo.hpp"

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
 Matrix * A,
 int rank,
 int iter_max, double res_tol, int monitor_iter,
 bool is_singular,
 bool diag_precon) 
  : Solver(), 
    A_(A),
    M_((diag_precon) ?
       (Matrix *)(new EnzoMatrixDiagonal) :
       (Matrix *)(new EnzoMatrixIdentity)),
    is_singular_(is_singular),
    rank_(rank),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    monitor_iter_(monitor_iter),
    rr0_(0),
    rr_min_(0),rr_max_(0),
    idensity_(0),  ipotential_(0),
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

  field_descr->ghost_depth    (idensity_,&gx_,&gy_,&gz_);

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  //  refresh(ir)->add_field(idensity_);
  refresh(ir)->add_all_fields(num_fields);

  /// Initialize matvec Refresh

#ifdef OLD_REFRESH
  id_refresh_matvec_ = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(id_refresh_matvec_)->add_all_fields(num_fields);
  //  refresh(id_refresh_matvec_)->add_field(ir_);
#endif

  TRACE_SOLVER("EnzoSolverCg() EXIT");
}

//----------------------------------------------------------------------

EnzoSolverCg::~EnzoSolverCg() throw ()
{
  INCOMPLETE("EnzoSolverCg::~EnzoSolverCg");
}

//----------------------------------------------------------------------

void EnzoSolverCg::pup (PUP::er &p)
{
  TRACEPUP;

  Solver::pup(p);

  p | A_;
  p | M_;
  p | is_singular_;
  p | rank_;
  p | iter_max_;
  p | res_tol_;
  p | monitor_iter_;
  p | rr0_;
  p | rr_min_;
  p | rr_max_;
  p | idensity_;
  p | ipotential_;
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
  INCOMPLETE("EnzoSolverCg::apply");
}
