// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverMg0.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoSolverMg0 class
///
/// Multigrid solver on a non-adaptive mesh.  Can be any mesh level, but
/// typically the root-grid (level = 0).
///
///======================================================================
///
///  "Coarse" view of MG0 multigrid solver
///
///   @code
///
///   $MG(A_h,X_h,B_h)$
///
///    while ( ! converged() ) 
///       if (level == min_level) then
///          solve_coarse()     solve $A_h X_h = B_h$
///       else
/// 1        p_pre_smooth()     smooth $A_h X_h = B_h$
/// 2        p_residual()       $R_h = B_h - A_h * X_h$
/// 3        p_restrict ()      $B_H = I_h^H R_h$
/// 4        MG()               solve $A_H X_H = B_H$  (repeat for W-cycle)
/// 5        p_prolong ()       $X_h = X_h + I_H^h X_H$
/// 6        p_post_smooth()    smooth $A_h X_h = B_h$
///
///  @endcode
///
///----------------------------------------------------------------------
///
///  "Fine" view of MG0 multigrid solver
///
///  @code
///
///  enter_solver()
///
///     iter = 0
///     initialize X,R,C
///     if (level == max_level) 
///        begin_cycle()
///
///  begin_cycle()
///
///     if (converged()) exit()
///     if (level == min_level) then
///        solve_coarse(A,X,B)
///     else
///        callback = p_pre_smooth()
///        call refresh (X,"level")
///
///  p_pre_smooth()
///
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
///
///  p_restrict_send(X)
///
///      A.residual(R,B,X) on level
///      pack R
///      index_parent.p_restrict_recv(R)
///
///  p_restrict_recv(B)
///
///      unpack B
///      --level
///      if (sync_restrict.next())
///          begin_cycle()
///      
///  coarse_solve(A,X,B)
///
///      solve A X = B
///      prolong_send(X)
///
///  prolong_send(X)
///
///      if (level < max_level)
///         for child           
///            pack X
///            child.prolong_recv(X)
///      else
///         begin_cycle()
///
///  prolong_recv(C)
///
///      ++level
///      unpack C         
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
///
///  p_post_smooth(A,X,B)
/// 
///      smooth.apply (A,X,B)
///      prolong_send()
///
///  @endcode
///
///======================================================================
///
/// Required Fields
///
/// - B                          linear system right-hand side
/// - R                          residual R = B - A*X
/// - X                          current solution to A*X = B
/// - C                          coarse-grid correction

#include "cello.hpp"
#include "charm_simulation.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#define NEW_REFRESH
// #define DEBUG_ENTRY
// #define DEBUG_SOLVER_MG0
// #define DEBUG_TRACE_LEVEL
// #define DEBUG_FIELD_MESSAGE
// #define DEBUG_SOLVER_INDEX

#ifdef DEBUG_TRACE_LEVEL
#   define TRACE_LEVEL(MSG,BLOCK)					\
  {									\
    int ia3[3];								\
    int it3[3];								\
    BLOCK->index().array(ia3,ia3+1,ia3+2);				\
    BLOCK->index().tree(it3,it3+1,it3+2);				\
    if (ia3[0]+ia3[1]+ia3[2]+it3[0]+it3[1]+it3[2]==0) {			\
      CkPrintf ("DEBUG_TRACE_LEVEL %d %s\n",BLOCK->level(),MSG);	\
    }									\
  }
#else
#   define TRACE_LEVEL(MSG,BLOCK) /*  this space for rent */
#endif  
  

#ifdef DEBUG_SOLVER_MG0
#   define TRACE_MG(block,msg)						\
  CkPrintf ("%d %s TRACE_MG %s\n",					\
	    CkMyPe(),(block != NULL) ? block->name().c_str() : "root",msg); \
  fflush(stdout);

#define DEBUG_FIELD(IX,NAME)						\
  {									\
    Data * data = enzo_block->data();					\
    Field field = data->field();					\
    T * X = (T*) field.values(IX);					\
    double xx=0.0;							\
    double yy=0.0;							\
    for (int iz=0; iz<mz_; iz++) {					\
      for (int iy=0; iy<my_; iy++) {					\
	for (int ix=0; ix<mx_; ix++) {					\
	  int i = ix + mx_*(iy + my_*iz);				\
	  xx+=X[i]*X[i];						\
	}								\
      }									\
    }									\
    for (int iz=gz_; iz<mz_-gz_; iz++) {				\
      for (int iy=gy_; iy<my_-gy_; iy++) {				\
	for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	  int i = ix + mx_*(iy + my_*iz);				\
	  yy+=X[i]*X[i];						\
	}								\
      }									\
    }									\
    CkPrintf ("%s:%d %s DEBUG_SOLVER ||%s|| = [%g] (%g)\n",		\
	      __FILE__,__LINE__,enzo_block->name().c_str(),NAME,xx,yy);	\
  }

#else
#   define TRACE_MG(block,msg) /*  ... */

#   define DEBUG_FIELD(IX,NAME) /* ... */
#endif

#ifdef DEBUG_FIELD_MESSAGE
#  define DEBUG_FIELD_MSG(BLOCK,MESSAGE)				\
  CkPrintf ("%d %s:%d %s DEBUG_FIELD_MESSAGE %s %d  %d %d %d  %lf  %lf\n", \
	    CkMyPe(),__FILE__,__LINE__,BLOCK->name().c_str(),		\
	    MESSAGE,							\
	    msg->n, msg->ic3[0],msg->ic3[1],msg->ic3[2],		\
	    msg->a[0],msg->a[msg->n-1]);					\
  fflush(stdout)
#else
#  define DEBUG_FIELD_MSG(BLOCK,MESSAGE) /* ...  */
#endif  

#ifdef DEBUG_SOLVER_INDEX
#  define DEBUG_INDEX(BLOCK,MESSAGE,INDX)			 \
  {									\
    int v3[3];								\
    BLOCK->index().values(v3);						\
    int w3[3] = {0};							\
    if ((INDX) != NULL) ((Index*)INDX)->values(w3);			\
    CkPrintf ("%d %s:%d  %08x %08x %08x --> %08x %08x %08x %s DEBUG_INDEX %s\n", \
	      CkMyPe(),__FILE__,__LINE__,				\
	      v3[0],v3[1],v3[2],					\
	      w3[0],w3[1],w3[2],					\
	      BLOCK->name().c_str(),					\
	      MESSAGE);							\
    fflush(stdout);							\
  }
#else
#  define DEBUG_INDEX(BLOCK,MESSAGE,INDEX)	\
  /* ... */
#endif
  
//======================================================================

EnzoSolverMg0::EnzoSolverMg0
(FieldDescr * field_descr, 
 int monitor_iter,
 int rank,
 int iter_max,
 int index_smooth_pre,
 int index_solve_coarse,
 int index_smooth_post,
 Restrict * restrict,
 Prolong * prolong,
 int min_level,
 int max_level) 
  : Solver(monitor_iter,min_level,max_level), 
    A_(NULL),
    index_smooth_pre_(index_smooth_pre),
    index_solve_coarse_(index_solve_coarse),
    index_smooth_post_(index_smooth_post),
    restrict_(restrict),
    prolong_(prolong),
    rank_(rank),
    iter_max_(iter_max), 
    ib_(0), ic_(0), ir_(0), ix_(0),
    mx_(0),my_(0),mz_(0),
    nx_(0),ny_(0),nz_(0),
    gx_(0),gy_(0),gz_(0),
    bs_(0), bc_(0)
    // [*]
{
  // Initialize temporary fields
  
  ir_ = field_descr->insert_temporary();
  ic_ = field_descr->insert_temporary();

  /// Initialize default Refresh

  add_refresh(4,0,neighbor_level,sync_barrier);

#ifdef NEW_REFRESH
  refresh(0)->add_all_fields();

  refresh(0)->add_field (ir_);
  refresh(0)->add_field (ic_);
#else
  refresh(0)->add_field (ir_);
  refresh(0)->add_field (ic_);
#endif  
}

//----------------------------------------------------------------------

EnzoSolverMg0::~EnzoSolverMg0 () throw()
// [*]
{
  delete prolong_;
  delete restrict_;

  prolong_ = NULL;
  restrict_ = NULL;
}


//----------------------------------------------------------------------

void EnzoSolverMg0::apply
( Matrix * A, int ix, int ib, Block * block) throw()
{

  TRACE_MG(block,"EnzoSolverMg0::apply()");

  Solver::begin_(block);

  A_ = A;
  ix_ = ix;
  ib_ = ib;

  Field field = block->data()->field();

  allocate_temporary_(field,block);
  
  TRACE_LEVEL("EnzoSolverMg0::apply",block);

  field.size           (&nx_,&ny_,&nz_);
  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  // Initialize child counter for restrict synchronization
  enzo_block->mg_sync_restrict_set_stop(NUM_CHILDREN(enzo_block->rank()));
  enzo_block->mg_sync_restrict_reset();
  enzo_block->mg_sync_prolong_set_stop(2); // self and parent
  enzo_block->mg_sync_prolong_reset();

  int precision = field.precision(ix_);

  if      (precision == precision_single)    
    enter_solver_<float>      (enzo_block);
  else if (precision == precision_double)    
    enter_solver_<double>     (enzo_block);
  else if (precision == precision_quadruple) 
    enter_solver_<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverMg0()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::enter_solver_ (EnzoBlock * enzo_block) throw()
/// [*]
///
///     iter = 0
///     initialize X,B,R,C
///     if (level == max_level) 
///        begin_cycle()
{
  TRACE_LEVEL("EnzoSolverMg0::enter_solver",enzo_block);

  TRACE_MG(enzo_block,"EnzoSolverMg0::enter_solver()");

  enzo_block->mg_iter_clear();

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X = (T*) field.values(ix_);
  T * R = (T*) field.values(ir_);
  T * C = (T*) field.values(ic_);

  // X = 0
  // R = B ( residual with X = 0 )
  // C = 0

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] = 0.0;
	R[i] = 0.0;
	C[i] = 0.0;
      }
    }
  }

  if (A_->is_singular()) {

    // Compute sum(B) and length() to project B onto range of A
    // if A is singular (

    long double reduce[2] = {0.0, 0.0};

    if (enzo_block->is_leaf()) {

      T* B = (T*) field.values(ib_);

      for (int iz=gz_; iz<nz_+gz_; iz++) {
	for (int iy=gy_; iy<ny_+gy_; iy++) {
	  for (int ix=gx_; ix<nx_+gx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    reduce[0] += B[i];
	  }
	}
      }
      
      reduce[1] = 1.0*nx_*ny_*nz_;
    }

#ifdef DEBUG_SOLVER_MG0    
    CkPrintf ("%d %s DEBUG_SOLVER_MG0 bs %lf bc %lf\n",
	      __LINE__,enzo_block->name().c_str(),double(reduce[0]),double(reduce[1]));
#endif    

    /// initiate callback for p_solver_mg0_shift_b and contribute to
    /// sum and count

    CkCallback callback(CkIndex_EnzoBlock::p_solver_mg0_shift_b(NULL), 
			enzo_block->proxy_array());

    enzo_block->contribute(2*sizeof(long double), &reduce, 
			   sum_long_double_2_type, callback);

  } else {
    begin_solve<T>(enzo_block);
  }

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_shift_b(CkReductionMsg* msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY enter p_solver_mg0_shift_b\n",
	      CkMyPe(),name().c_str(),this);
#endif

  /// EnzoBlock accumulates global contributions to SUM(B) and COUNT(B)
  EnzoSolverMg0* solver = 
    static_cast<EnzoSolverMg0*> (this->solver());
  
  long double* data = (long double*) msg->getData();

  solver->set_bs( data[0] );
  solver->set_bc( data[1] );

#ifdef DEBUG_SOLVER_MG0  
  CkPrintf ("%d %s DEBUG_SOLVER_MG0 bs %lf bc %lf\n",
	    __LINE__,name().c_str(),double(data[0]),double(data[1]));
#endif  

  delete msg;

  /// Start solver
  Field field = this->data()->field();
  int precision = field.precision(0);
  if      (precision == precision_single)    
    solver->begin_solve<float>(this);
  else if      (precision == precision_double)    
    solver->begin_solve<double>(this);
  else if      (precision == precision_quadruple)    
    solver->begin_solve<long double>(this);
  else 
    ERROR1("EnzoSolverMg0::p_solver_mg0_shift_b()",
	   "precision %d not recognized", precision);

#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY  exit p_solver_mg0_shift_b\n",
	      CkMyPe(),name().c_str(),this);
#endif
    performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::begin_solve(EnzoBlock * enzo_block) throw()
{
  TRACE_LEVEL("EnzoSolverMg0::begin_solve",enzo_block);
  // start the MG V-cycle with the max level blocks


  if (A_->is_singular() && (enzo_block->level() == max_level_)) {

    // Shift B if needed to be in range(A) for periodic b.c.
    
    Field field = enzo_block->data()->field();
    T* B  = (T*) field.values(ib_);
    T shift = -bs_ / bc_;
#ifdef DEBUG_SOLVER_MG0      
    CkPrintf ("%d %s DEBUG_SOLVER_MG0 bs %lf bc %lf\n",
	      __LINE__,enzo_block->name().c_str(),double(bs_),double(bc_));
#endif      

    TRACE_FIELD_("B shift 0",B,1.0);
    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  B[i] += shift;
	}
      }
    }

    TRACE_FIELD_("B shift 1",B,1.0);
  }

  // control flow starts at leaves, even in level > max_level,
  // since coarse solve may require reductions over all Blocks

  if (enzo_block->is_leaf()) {
    begin_cycle_<T>(enzo_block);
  }
  
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::begin_cycle_(EnzoBlock * enzo_block) throw()
/// [*]
///
///     if (converged()) exit()
///     if (level == min_level) then
///        coarse_solve(A,X,B)
///     else
///        callback = p_pre_smooth()
///        call refresh (X,"level")
{
  TRACE_LEVEL("EnzoSolverMg0::begin_cycle",enzo_block);
  const int level = enzo_block->level();

  TRACE_MG(enzo_block,"EnzoSolverMg0::begin_cycle()");

  // Monitor output
  
  bool is_converged = is_converged_(enzo_block);

  const int iter = enzo_block->mg_iter();

  const bool l_output =
    ( ( enzo_block->index().is_zero() && level == max_level_) &&
      ( (iter == 0) ||
	(is_converged) ||
	(monitor_iter_ && (iter % monitor_iter_) == 0 )) );

  if (l_output) {
    Monitor* monitor = enzo_block->simulation()->monitor();
    monitor_output_(enzo_block,iter,0.0,0.0,0.0,0.0);
  }

  if (is_converged) {

    TRACE_MG(enzo_block,"converged");

    // CkCallback(callback_,
    // 	       CkArrayIndexIndex(enzo_block->index()),
    // 	       enzo_block->proxy_array()).send();

  } else if (level == min_level_) {
    // TODO REFRESH X
    TRACE_MG(enzo_block,"calling coarse solve");

    Field field = enzo_block->data()->field();
    T * X = (T*) field.values(ix_);

    std::fill_n(X,mx_*my_*mz_,0.0);

    Simulation * simulation = proxy_simulation.ckLocalBranch();
    Solver * solve_coarse = simulation->problem()->solver(index_solve_coarse_);

    solve_coarse->set_min_level(min_level_);
    solve_coarse->set_max_level(min_level_);
    solve_coarse->set_sync_id (4);
    solve_coarse->set_callback(CkIndex_EnzoBlock::p_solver_mg0_solve_coarse());
  
    solve_coarse->apply(A_,ix_,ib_,enzo_block);

  } else {

    TRACE_MG(enzo_block,"calling smoother");

    if ( (! enzo_block->is_leaf()) && (level < max_level_) ) {

      Field field = enzo_block->data()->field();
      T * X = (T*) field.values(ix_);
      std::fill_n(X,mx_*my_*mz_,0.0);

    }

#ifdef DEBUG_SOLVER_REFRESH    
    CkPrintf ("DEBUG_SOLVER_MG refresh sync_face %d\n",refresh.sync_id());
#endif

    if (index_smooth_pre_ >= 0) {

      Data * data = enzo_block->data();
      Field field = data->field();

      Simulation * simulation = proxy_simulation.ckLocalBranch();

      Solver * smooth_pre = simulation->problem()->solver(index_smooth_pre_);

      smooth_pre->set_min_level(enzo_block->level());
      smooth_pre->set_max_level(enzo_block->level());
      smooth_pre->set_sync_id (6);
      smooth_pre->set_callback(CkIndex_EnzoBlock::p_solver_mg0_pre_smooth());

      smooth_pre->apply(A_,ix_,ib_,enzo_block);

    } else {

      pre_smooth<T>(enzo_block);

    }

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_solve_coarse()
/// [*]
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  long double reduce[1] = {0};
  
  CkCallback callback(CkIndex_EnzoBlock::p_solver_mg0_barrier(NULL), 
		      proxy_array());
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY before barrier\n",
	      CkMyPe(),name().c_str(),this);
#endif
  contribute(0, reduce, sum_long_double_type, callback);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_barrier(CkReductionMsg* msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  delete msg;

#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY  after barrier\n",
	      CkMyPe(),name().c_str(),this);
#endif
  TRACE_MG(this,"EnzoBlock::solver_mg0_coarse_solve()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  
  // assume all fields have same precision
  int precision = field.precision(0);
  if      (precision == precision_single)    
    solver->solve_coarse<float>(enzo_block);
  else if (precision == precision_double)    
    solver->solve_coarse<double>(enzo_block);
  else if (precision == precision_quadruple) 
    solver->solve_coarse<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverMg0::p_solver_mg0_solve_coarse()",
	   "precision %d not recognized", precision);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_pre_smooth()
/// [*]
{
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY enter p_solver_mg0_pre_smooth\n",
	      CkMyPe(),name().c_str(),this);
#endif
  TRACE_MG(this,"EnzoBlock::solver_mg0_pre_smooth()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();

  // assume all fields have same precision

  int precision = field.precision(0);
  if      (precision == precision_single)    
    solver->pre_smooth<float>(enzo_block);
  else if (precision == precision_double)    
    solver->pre_smooth<double>(enzo_block);
  else if (precision == precision_quadruple) 
    solver->pre_smooth<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverMg0::p_solver_mg0_pre_smooth()",
	   "precision %d not recognized", precision);
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY  exit p_solver_mg0_pre_smooth\n",
	      CkMyPe(),name().c_str(),this);
#endif
    performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::pre_smooth(EnzoBlock * enzo_block) throw()
/// [*]
///
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  TRACE_LEVEL("EnzoSolverMg0::pre_smooth",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::pre_smooth()");

  restrict_send<T>(enzo_block);

  // All Blocks must call coarse solver since may involve
  // global reductions
  
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  Solver * solve_coarse = simulation->problem()->solver(index_solve_coarse_);

  solve_coarse->set_min_level(min_level_);
  solve_coarse->set_max_level(min_level_);
  solve_coarse->set_sync_id (4);
  solve_coarse->set_callback(CkIndex_EnzoBlock::p_solver_mg0_solve_coarse());
  
  //  TRACE_FIELD_("X",X,1.0);
  //  TRACE_FIELD_("B",B,1.0);
  solve_coarse->apply(A_,ix_,ib_,enzo_block);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::restrict_send(EnzoBlock * enzo_block) throw()
/// 
/// [*] restrict send
///
///      A.residual(R,B,X)
///      pack R
///      index_parent.p_restrict_recv(R)
{
  TRACE_LEVEL("EnzoSolverMg0::restrict_send",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_send()");

  Data * data = enzo_block->data();
  Field field = data->field();

  A_->residual(ir_, ib_, ix_, enzo_block);

  Index index        = enzo_block->index();
  Index index_parent = index.index_parent(min_level_);
  const  int level   = index.level();  
  // copy face data to FieldFace

  // Pack and send "R" to parent

  int ic3[3];
  index.child(level,&ic3[0],&ic3[1],&ic3[2],min_level_);
  
  // <COMMON CODE> in restrict_send_() and prolong_send_()
  
  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  Refresh * refresh = new Refresh;
  refresh->add_field(ir_);

  // copy data from EnzoBlock to array via FieldFace

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, lg3, refresh_coarse, refresh, true);

  field_face->set_restrict(restrict_);
  
  int narray; 
  char * array;

  field_face->face_to_array(enzo_block->data()->field(),&narray,&array);

  delete field_face;

  // Create a FieldMsg for sending data to parent
  // (note: charm messages not deleted on send; are deleted on receive)

  FieldMsg * msg  = new (narray) FieldMsg;
 
  /// WARNING: double copy

  // Copy FieldFace data to msg

  msg->n = narray;
  memcpy (msg->a, array, narray);
  delete [] array;
  msg->ic3[0] = ic3[0];
  msg->ic3[1] = ic3[1];
  msg->ic3[2] = ic3[2];

  //  </COMMON CODE>

  //  TRACE_FIELD_("B",B,1.0);
  //  TRACE_FIELD_("X",X,1.0); // XXXXX
  //  TRACE_FIELD_("R",R,1.0); // XXXXX

  DEBUG_FIELD_MSG(enzo_block,"restrict_send");
  DEBUG_INDEX(enzo_block,"restrict_send",(&index_parent));
  
  enzo_block->thisProxy[index_parent].p_solver_mg0_restrict_recv(msg);

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_restrict_recv(FieldMsg * msg)
/// [*]
{
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY enter p_solver_mg0_restrict_recv\n",
	      CkMyPe(),name().c_str(),this);
#endif
  DEBUG_FIELD_MSG(this,"restrict_inter");

  TRACE_MG(this,"EnzoBlock::restrict_recv()");
  DEBUG_INDEX(this,"restrict_inter",NULL);

  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  Field field = data()->field();
  // assume all fields have same precision
  int precision = field.precision(0);
  if      (precision == precision_single)    
    solver->restrict_recv<float>(this,msg);
  else if (precision == precision_double)    
    solver->restrict_recv<double>(this,msg);
  else if (precision == precision_quadruple) 
    solver->restrict_recv<long double>(this,msg);
  else 
    ERROR1("EnzoSolverMg0::p_solver_mg0_restrict_recv()",
	   "precision %d not recognized", precision);
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY  exit p_solver_mg0_restrict_recv\n",
	      CkMyPe(),name().c_str(),this);
#endif
    performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::restrict_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
/// 
/// [*] restrict recv
///
///      [ unpack B ]
///      if (sync.next())
///          begin_cycle()
{

  DEBUG_FIELD_MSG(enzo_block,"restrict_recv");

  // Unpack "B" vector data from children

  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  Refresh * refresh = new Refresh;
  refresh->add_field(ib_);

  // copy data from msg to this EnzoBlock

  int * ic3 = msg->ic3;

  FieldFace * field_face = enzo_block->create_face 
    (if3, ic3, lg3, refresh_coarse, refresh, true);

  field_face->set_restrict(restrict());

  char * a = msg->a;
  field_face->array_to_face(a, enzo_block->data()->field());
  delete field_face;

  delete msg;

  // continue with EnzoSolverMg0

  TRACE_LEVEL("EnzoSolverMg0::restrict_recv",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_recv()");

  //  TRACE_FIELD_("B",B,1.0);
  
  if (enzo_block->mg_sync_restrict_next()) {
    begin_cycle_<T>(enzo_block);
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::solve_coarse(EnzoBlock * enzo_block) throw()
/// 
/// [*] coarse solve
///
///      solve A X = B
///      end_cycle()
{
  TRACE_LEVEL("EnzoSolverMg0::solve_coarse",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::solver_coarse()");
 
  /// Prolong solution to next-finer level

  const int level = enzo_block->level();
  
  if (level == min_level_) {

    if ( (! enzo_block->is_leaf()) && (level < max_level_) ) {

      prolong_send_<T>(enzo_block);
      
    }

    end_cycle<T>(enzo_block);

  } else if (level > min_level_) {

    enzo_block->p_solver_mg0_prolong_recv(NULL);

  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::prolong_send_(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong send
///
///      for child           
///         pack X
///         child.prolong_recv(X)
{
  TRACE_LEVEL("EnzoSolverMg0::prolong_send",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::prolong_send()");

  ItChild it_child(enzo_block->rank());
  int ic3[3];
  //  TRACE_FIELD_("X",X,1.0);
  while (it_child.next(ic3)) {

    Index index_child = enzo_block->index().index_child(ic3,min_level_);

    // Pack and send "X" to children

    // <COMMON CODE> in restrict_send_() and prolong_send_()

    int if3[3] = {0,0,0};
    bool lg3[3] = {false,false,false};
    Refresh * refresh = new Refresh;
    refresh->add_field(ix_);
    
    // copy data from EnzoBlock to array via FieldFace

    FieldFace * field_face = enzo_block->create_face
      (if3, ic3, lg3, refresh_fine, refresh, true);

    field_face->set_prolong(prolong_);

    int narray; 
    char * array;
    
    field_face->face_to_array (enzo_block->data()->field(),&narray,&array);

    delete field_face;

    // Create a FieldMsg for sending data to parent
    // (note: charm messages not deleted on send; are deleted on receive)
    
    FieldMsg * msg  = new (narray) FieldMsg;

    /// WARNING: double copy

    // Copy FieldFace data to msg

    msg->n = narray;
    memcpy (msg->a, array, narray);
    delete [] array;
    msg->ic3[0] = ic3[0];
    msg->ic3[1] = ic3[1];
    msg->ic3[2] = ic3[2];

    //  </COMMON CODE>

    DEBUG_FIELD_MSG(enzo_block,"prolong_send");

    DEBUG_INDEX(enzo_block,"prolong_send",&index_child);

    enzo_block->thisProxy[index_child].p_solver_mg0_prolong_recv(msg);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  // Save message
  if (msg != NULL) mg_msg_ = msg;

  // Return if not ready yet
  if (! mg_sync_prolong_next())
    return;

  // Restore saved message
  msg = mg_msg_;
  mg_msg_ = NULL;
  
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY enter p_solver_mg0_prolong_recv\n",
	      CkMyPe(),name().c_str(),this);
#endif
  DEBUG_INDEX(this,"prolong_inter",NULL);

  TRACE_MG (this,"EnzoBlock::p_solver_mg0_prolong_recv()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  Field field = data()->field();
  // assume all fields have same precision
  int precision = field.precision(0);
  if      (precision == precision_single)    
    solver->prolong_recv<float>(this,msg);
  else if (precision == precision_double)    
    solver->prolong_recv<double>(this,msg);
  else if (precision == precision_quadruple) 
    solver->prolong_recv<long double>(this,msg);
  else 
    ERROR1("EnzoSolverMg0::p_solver_mg0_prolong_recv()",
	   "precision %d not recognized", precision);
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY  exit p_solver_mg0_prolong_recv\n",
	      CkMyPe(),name().c_str(),this);
#endif
    performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::prolong_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
/// 
/// [ ] prolong recv
///
///      [ unpack C ]
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{

  // Unpack "C" vector data from children

  DEBUG_FIELD_MSG(enzo_block,"prolong_recv");
    
  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  Refresh * refresh = new Refresh;
  refresh->add_field(ic_);

  // copy data from msg to this EnzoBlock

  FieldFace * field_face = enzo_block->create_face 
    (if3, msg->ic3, lg3, refresh_fine, refresh, true);

  field_face->set_prolong(prolong());

  field_face->array_to_face (msg->a, enzo_block->data()->field());

  delete field_face;

  // GET DATA
  
  delete msg;

  TRACE_LEVEL("EnzoSolverMg0::prolong_recv",enzo_block);
  TRACE_MG (enzo_block,"EnzoSolverMg0::prolong_recv()");
  
  Data * data = enzo_block->data();
  Field field = data->field();

  T * X = (T*) field.values(ix_);
  T * C = (T*) field.values(ic_);

  TRACE_FIELD_("C",C,1.0);

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] += C[i];
      }
    }
  }

  if (index_smooth_post_ >= 0) {
    Simulation * simulation = proxy_simulation.ckLocalBranch();
    Solver * smooth_post = simulation->problem()->solver(index_smooth_post_);

    smooth_post->set_min_level(enzo_block->level());
    smooth_post->set_max_level(enzo_block->level());
    smooth_post->set_sync_id (8);
    smooth_post->set_callback(CkIndex_EnzoBlock::p_solver_mg0_post_smooth());
  
    smooth_post->apply(A_,ix_,ib_,enzo_block);

  } else {

    post_smooth<T>(enzo_block);
  }

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_post_smooth()
/// [*]
{
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY enter p_solver_mg0_post_smooth\n",
	      CkMyPe(),name().c_str(),this);
#endif
  TRACE_MG (this,"EnzoBlock::p_solver_mg0_post_smooth()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  
  // assume all fields have same precision
  int precision = field.precision(0);
  if      (precision == precision_single)    
    solver->post_smooth<float>(enzo_block);
  else if (precision == precision_double)    
    solver->post_smooth<double>(enzo_block);
  else if (precision == precision_quadruple) 
    solver->post_smooth<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverMg0::p_solver_mg0_post_smooth()",
	   "precision %d not recognized", precision);

#ifdef DEBUG_ENTRY
    CkPrintf ("%d %s %p mg0 DEBUG_ENTRY  exit p_solver_mg0_post_smooth\n",
	      CkMyPe(),name().c_str(),this);
#endif
    performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::post_smooth(EnzoBlock * enzo_block) throw()
/// post smooth
/// [*]
///
///      smooth.apply (A,X,B)
///      end_cycle()
{
  TRACE_LEVEL("EnzoSolverMg0::post_smooth",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::post_smooth()");

  const int level = enzo_block->level();

  if ( ( ! enzo_block->is_leaf() ) && (level < max_level_)) {

    prolong_send_<T>(enzo_block);
  } 

  TRACE_MG(enzo_block,"EnzoSolverMg0::post_smooth() calling end_cycle()");
  end_cycle<T>(enzo_block);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::end_cycle(EnzoBlock * enzo_block) throw()
/// 
/// [X] end cycle
///
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
{
  TRACE_LEVEL("EnzoSolverMg0::end_cycle",enzo_block);

  TRACE_MG(enzo_block,"EnzoSolverMg0::end_cycle()");
  
  enzo_block->mg_iter_increment();

  bool is_converged = is_converged_(enzo_block);

  const int iter = enzo_block->mg_iter();
	    
  const int level = enzo_block->level();

  const bool l_output =
    ( ( enzo_block->index().is_zero() && level == max_level_) &&
      ( (is_converged) ||
	(monitor_iter_ && (iter % monitor_iter_) == 0 )) );

  if (l_output) {
    Monitor* monitor = enzo_block->simulation()->monitor();
    monitor_output_(enzo_block,iter,0.0,0.0,0.0,0.0);
  }

  if (is_converged) {

    end_(enzo_block);

  } else if (enzo_block->is_leaf() || (enzo_block->level() == max_level_)) {

    begin_cycle_<T>(enzo_block);

  }
}


//======================================================================

bool EnzoSolverMg0::is_converged_(EnzoBlock * enzo_block) const
/// [*]
{
  TRACE_MG(enzo_block,"EnzoSolverMg0::is_converged");
  return (enzo_block->mg_iter() >= iter_max_);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::end_(Block * block)
{
  Field field = block->data()->field();

  deallocate_temporary_(field,block);
    
  Solver::end_(block);

  CkCallback(callback_,
	     CkArrayIndexIndex(block->index()),
	     block->proxy_array()).send();

}
