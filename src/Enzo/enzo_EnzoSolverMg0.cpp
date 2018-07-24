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

// #define DEBUG_COPY

// #define DEBUG_NEGATE_X
// #define DEBUG_ENTRY
// #define DEBUG_SOLVER_MG0
// #define DEBUG_SOLVER_CONTROL
// #define DEBUG_TRACE_LEVEL
// #define DEBUG_FIELD_MESSAGE
// #define DEBUG_SOLVER_INDEX
// #define DEBUG_PROLONG

#define AFTER_CYCLE(BLOCK,CYCLE) (BLOCK->cycle() >= CYCLE)
#define CYCLE 0

#ifdef DEBUG_TRACE_LEVEL
#   define TRACE_LEVEL(MSG,BLOCK)					\
  {									\
    if (AFTER_CYCLE(BLOCK,CYCLE))					\
      CkPrintf ("%s %s DEBUG_TRACE_LEVEL %d %s\n",BLOCK->name().c_str(), \
		name_.c_str(),BLOCK->level(),MSG);			\
  }
#else
#   define TRACE_LEVEL(MSG,BLOCK) /*  this space for rent */
#endif  
  

#ifdef DEBUG_SOLVER_MG0
#   define TRACE_MG(block,msg)						\
  if (AFTER_CYCLE(block,CYCLE))						\
    CkPrintf ("%d %s TRACE_MG %s %s\n",					\
	      CkMyPe(),(block != NULL) ? block->name().c_str() : "root",name_.c_str(),msg); \
  fflush(stdout);
#   define TRACE_MG_BLOCK(block,msg)					\
  if (AFTER_CYCLE(block,CYCLE))						\
    CkPrintf ("%d %s TRACE_MG %s\n",					\
	      CkMyPe(),(block != NULL) ? block->name().c_str() : "root",msg); \
  fflush(stdout);

#define DEBUG_FIELD(BLOCK,IX,NAME)			\
  {							\
    if (AFTER_CYCLE(BLOCK,CYCLE)) {			\
      Data * data = enzo_block->data();			\
      Field field = data->field();			\
      enzo_float * X = (enzo_float*) field.values(IX);	\
      double xx=0.0;					\
      double yy=0.0;					\
      for (int iz=0; iz<mz_; iz++) {			\
	for (int iy=0; iy<my_; iy++) {			\
	  for (int ix=0; ix<mx_; ix++) {		\
	    int i = ix + mx_*(iy + my_*iz);		\
	    xx+=X[i]*X[i];				\
	  }						\
	}						\
      }							\
      for (int iz=gz_; iz<mz_-gz_; iz++) {		\
	for (int iy=gy_; iy<my_-gy_; iy++) {		\
	  for (int ix=gx_; ix<mx_-gx_; ix++) {		\
	    int i = ix + mx_*(iy + my_*iz);		\
	    yy+=X[i]*X[i];				\
	  }						\
	}						\
      }							\
      CkPrintf ("%s DEBUG_FIELD ||%s|| = [%g] (%g)\n",	\
		enzo_block->name().c_str(),NAME,xx,yy);	\
    }							\
  }

#else
#   define TRACE_MG(block,msg) /*  ... */
#   define TRACE_MG_BLOCK(block,msg) /*  ... */

#   define DEBUG_FIELD(BLOCK,IX,NAME) /* ... */
#endif

#ifdef DEBUG_FIELD_MESSAGE
#  define DEBUG_FIELD_MSG(BLOCK,MESSAGE)				\
  if (AFTER_CYCLE(BLOCK,CYCLE))						\
    CkPrintf ("%d %s %s DEBUG_FIELD_MESSAGE %s %d  %d %d %d  %lf  %lf\n", \
	      CkMyPe(),__FILE__,BLOCK->name().c_str(),			\
	      MESSAGE,							\
	      msg->n, msg->ic3[0],msg->ic3[1],msg->ic3[2],		\
	      msg->a[0],msg->a[msg->n-1]);				\
  fflush(stdout)
#else
#  define DEBUG_FIELD_MSG(BLOCK,MESSAGE) /* ...  */
#endif  

#ifdef DEBUG_SOLVER_INDEX
#  define DEBUG_INDEX(BLOCK,MESSAGE,INDX)				\
  {									\
    int v3[3];								\
    BLOCK->index().values(v3);						\
    int w3[3] = {0};							\
    if ((INDX) != NULL) ((Index*)INDX)->values(w3);			\
    if (AFTER_CYCLE(BLOCK,CYCLE))					\
      CkPrintf ("%d %s  %08x %08x %08x --> %08x %08x %08x %s DEBUG_INDEX %s\n", \
		CkMyPe(),__FILE__,					\
		v3[0],v3[1],v3[2],					\
		w3[0],w3[1],w3[2],					\
		BLOCK->name().c_str(),					\
		MESSAGE);						\
    fflush(stdout);							\
  }
#else
#  define DEBUG_INDEX(BLOCK,MESSAGE,INDEX)	\
  /* ... */
#endif


#ifdef DEBUG_SOLVER_CONTROL
#   define SOLVER_CONTROL(BLOCK,MIN,MAX,MESSAGE)			\
  if (AFTER_CYCLE(BLOCK,CYCLE))						\
    CkPrintf ("DEBUG_SOLVER_CONTROL %4s %8s <= L <= %8s %30s %20s\n",	\
	      name_.c_str(),MIN,MAX,MESSAGE,BLOCK->name().c_str());
#else
#   define SOLVER_CONTROL(BLOCK,MIN,MAX,MESSAGE) /* ... */
#endif
  
//======================================================================

EnzoSolverMg0::EnzoSolverMg0
(std::string name,
 int monitor_iter,
 int restart_cycle,
 int rank,
 int iter_max,
 double res_tol,
 int index_smooth_pre,
 int index_solve_coarse,
 int index_smooth_post,
 int index_smooth_last,
 Restrict * restrict,
 Prolong * prolong,
 int min_level,
 int max_level,
 int min_level_coarse,
 int max_level_coarse,
 bool is_unigrid) 
  : Solver(name,monitor_iter,restart_cycle,min_level,max_level,is_unigrid),
    min_level_coarse_(min_level_coarse),
    max_level_coarse_(max_level_coarse),
    A_(NULL),
    index_smooth_pre_(index_smooth_pre),
    index_solve_coarse_(index_solve_coarse),
    index_smooth_post_(index_smooth_post),
    index_smooth_last_(index_smooth_last),
    restrict_(restrict),
    prolong_(prolong),
    rank_(rank),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    ib_(0), ic_(0), ir_(0), ix_(0),
    mx_(0),my_(0),mz_(0),
    nx_(0),ny_(0),nz_(0),
    gx_(0),gy_(0),gz_(0),
    bs_(0), bc_(0),
    rr_(0), rr_local_(0), rr0_(0)
{
  // Initialize temporary fields

  FieldDescr * field_descr = cello::field_descr();
  
  ir_ = field_descr->insert_temporary();
  ic_ = field_descr->insert_temporary();

  /// Initialize default Refresh

  add_refresh(4,0,neighbor_level,sync_barrier,
	      enzo_sync_id_solver_mg0);

  refresh(0)->add_field (ir_);
  refresh(0)->add_field (ic_);

#ifdef NEW_SYNC  
  ScalarDescr * scalar_descr_int  = cello::scalar_descr_int();
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();

  iscalar_iter_   = scalar_descr_int->new_value(name + ":iter");
  isync_restrict_ = scalar_descr_sync->new_value(name + ":restrict");
  isync_prolong_  = scalar_descr_sync->new_value(name + ":prolong");
#endif  

}

//----------------------------------------------------------------------

EnzoSolverMg0::~EnzoSolverMg0 () throw()
{
  delete prolong_;
  delete restrict_;

  prolong_ = NULL;
  restrict_ = NULL;
}


//----------------------------------------------------------------------

void EnzoSolverMg0::apply
( std::shared_ptr<Matrix> A, int ix, int ib, Block * block) throw()
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


#ifdef NEW_SYNC  
  ScalarData<Sync> * scalar_data  = block->data()->scalar_data_sync();
  ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
  
  Sync & sync_restrict = scalar_data->value(scalar_descr,isync_restrict_);

  sync_restrict.set_stop(NUM_CHILDREN(enzo_block->rank()));
  sync_restrict.reset();
  
  Sync & sync_prolong = scalar_data->value(scalar_descr,isync_prolong_);
  sync_prolong.set_stop(2); // self and parent
  sync_prolong.reset();

#else  
  
  enzo_block->mg_sync_restrict_set_stop(NUM_CHILDREN(enzo_block->rank()));
  enzo_block->mg_sync_restrict_reset();
  enzo_block->mg_sync_prolong_set_stop(2); // self and parent
  enzo_block->mg_sync_prolong_reset();
#endif
  
#ifdef DEBUG_PROLONG
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%s %s %s DEBUG_PROLONG prolong_set_stop(2)\n",
	      __FILE__,block->name().c_str(),name_.c_str());
#endif  

  enter_solver_ (enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::enter_solver_ (EnzoBlock * enzo_block) throw()
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

  enzo_float * X = (enzo_float*) field.values(ix_);
  enzo_float * R = (enzo_float*) field.values(ir_);
  enzo_float * C = (enzo_float*) field.values(ic_);

  // X = 0
  // R = B ( residual with X = 0 )
  // C = 0

  std::fill_n(X,mx_*my_*mz_,0.0);
  std::fill_n(R,mx_*my_*mz_,0.0);
  std::fill_n(C,mx_*my_*mz_,0.0);

  if (A_->is_singular()) {

    // Compute sum(B) and length() to project B onto range of A
    // if A is singular (

    long double reduce[2] = {0.0, 0.0};

    if (is_finest_(enzo_block)) {

      enzo_float* B = (enzo_float*) field.values(ib_);

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
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%s DEBUG_SOLVER_MG0 bs %lf bc %lf\n",
		enzo_block->name().c_str(),double(reduce[0]),double(reduce[1]));
#endif    

    /// initiate callback for p_solver_mg0_shift_b and contribute to
    /// sum and count

    DEBUG_FIELD(enzo_block,ix_,"X");
    DEBUG_FIELD(enzo_block,ib_,"B");
    CkCallback callback(CkIndex_EnzoBlock::r_solver_mg0_shift_b(NULL), 
			enzo_block->proxy_array());

    SOLVER_CONTROL (enzo_block,"min","max","1A calling begin_solve (shift)");

    enzo_block->contribute(2*sizeof(long double), &reduce, 
			   sum_long_double_2_type, callback);

  } else {

    DEBUG_FIELD(enzo_block,ix_,"X");
    DEBUG_FIELD(enzo_block,ib_,"B");

    SOLVER_CONTROL(enzo_block,"min","max","1B calling begin_solve (no shift)");

    begin_solve (enzo_block,NULL);

  }

}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_mg0_shift_b(CkReductionMsg* msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverMg0*> (solver())->begin_solve(this,msg);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::begin_solve(EnzoBlock * enzo_block,
				CkReductionMsg *msg) throw()
{
  TRACE_LEVEL("EnzoSolverMg0::begin_solve",enzo_block);
  // start the MG V-cycle with the max level blocks


  SOLVER_CONTROL(enzo_block,"min","max", "2 setting shift");
  
  if (msg != NULL) {
    
    long double* data = (long double*) msg->getData();

    bs_ = data[0];
    bc_ = data[1];

    delete msg;
  } 
  
  if (A_->is_singular() && is_finest_(enzo_block)) {

    // Shift B if needed to be in range(A) for periodic b.c.
    
    SOLVER_CONTROL(enzo_block,"fine","fine", "3 shifting B");

    Field field = enzo_block->data()->field();
    enzo_float* B  = (enzo_float*) field.values(ib_);
    enzo_float shift = -bs_ / bc_;
#ifdef DEBUG_SOLVER_MG0      
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%s DEBUG_SOLVER_MG0 bs B %lf bc %lf\n",
		enzo_block->name().c_str(),double(bs_),double(bc_));
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

    DEBUG_FIELD (enzo_block,ib_,"B shift");
  }

  // control flow starts at leaves, even in level > max_level,
  // since coarse solve may require reductions over all Blocks

  if (is_finest_(enzo_block)) {

    SOLVER_CONTROL(enzo_block, "fine","max", "4 begin V-cycle");

    begin_cycle_ (enzo_block);
  }
  
}

//----------------------------------------------------------------------

void EnzoSolverMg0::begin_cycle_(EnzoBlock * enzo_block) throw()
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
  
  // bool is_converged = is_converged_(enzo_block);

  const int iter = enzo_block->mg_iter();

  const bool l_output =
    ( ( enzo_block->index().is_zero() && level == max_level_) &&
      ( (iter == 0))); // ||

  if (l_output) {
    monitor_output_(enzo_block,iter,rr0_,0.0,rr_,0.0);
  }

  // if (is_converged) {
  
  TRACE_MG(enzo_block,"converged");

  // CkCallback(callback_,
  // 	       CkArrayIndexIndex(enzo_block->index()),
  // 	       enzo_block->proxy_array()).send();

  if (level <= min_level_) {

    SOLVER_CONTROL(enzo_block,"min","coarse","5B coarse solve V-cycle");

    //    if (!is_active_(enzo_block)) {
    // TODO REFRESH X
    TRACE_MG(enzo_block,"calling coarse solve");

    Field field = enzo_block->data()->field();
    enzo_float * X = (enzo_float*) field.values(ix_);

    std::fill_n(X,mx_*my_*mz_,0.0);

    call_coarse_solver(enzo_block);

  } else {

    SOLVER_CONTROL(enzo_block,"coarse+1","fine", "5A pre-smooth");

    TRACE_MG(enzo_block,"calling smoother");

    if ( ! is_finest_(enzo_block) ) {

      Field field = enzo_block->data()->field();
      enzo_float * X = (enzo_float*) field.values(ix_);
      std::fill_n(X,mx_*my_*mz_,0.0);

    }

#ifdef DEBUG_SOLVER_REFRESH    
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("DEBUG_SOLVER_MG refresh sync_face %d\n",refresh.sync_id());
#endif

    if (index_smooth_pre_ >= 0) {

      Data * data = enzo_block->data();
      Field field = data->field();

      Simulation * simulation = proxy_simulation.ckLocalBranch();

      Solver * smooth_pre = simulation->problem()->solver(index_smooth_pre_);

      smooth_pre->set_min_level(enzo_block->level());
      smooth_pre->set_max_level(enzo_block->level());
      smooth_pre->set_sync_id (enzo_sync_id_solver_mg0_pre);
      smooth_pre->set_callback(CkIndex_EnzoBlock::p_solver_mg0_pre_smooth());

      smooth_pre->apply(A_,ix_,ib_,enzo_block);

    } else {

      pre_smooth (enzo_block);

    }

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_solve_coarse()
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY before barrier\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

#ifdef DEBUG_SOLVER_MG0
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("DEBUG_SOLVER_MG0 %s solver->rr_local() = %llg\n",
	      name().c_str(),solver->rr_local());
#endif

  CkCallback callback(CkIndex_EnzoBlock::r_solver_mg0_barrier(NULL), 
		      proxy_array());
  long double data[1] = {solver->rr_local()};
  contribute(sizeof(long double), data, sum_long_double_type, callback);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_mg0_barrier(CkReductionMsg* msg)
{
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

#ifdef DEBUG_SOLVER_MG0
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("DEBUG_SOLVER_MG0 %s rr_ %llg\n",
	      name().c_str(),solver->rr());
#endif

  performance_start_(perf_compute,__FILE__,__LINE__);

  long double rr = ((long double*) msg->getData())[0];
  solver->set_rr(rr);
  solver->set_rr_local(0.0);
  if (mg_iter_==0) solver->set_rr0(rr);
  
  delete msg;

#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  after barrier\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  TRACE_MG_BLOCK(this,"EnzoBlock::solver_mg0_coarse_solve()");
  
  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  solver->solve_coarse(enzo_block);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_pre_smooth()
{
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY enter p_solver_mg0_pre_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  TRACE_MG_BLOCK(this,"EnzoBlock::solver_mg0_pre_smooth()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  solver->pre_smooth(enzo_block);

#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  exit p_solver_mg0_pre_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::pre_smooth(EnzoBlock * enzo_block) throw()
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "6 pre-smooth");
  
  TRACE_LEVEL("EnzoSolverMg0::pre_smooth",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::pre_smooth()");

  restrict_send (enzo_block);

  // All Blocks must call coarse solver since may involve
  // global reductions

  call_coarse_solver(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_coarse_solver(EnzoBlock * enzo_block) throw()
{
  SOLVER_CONTROL(enzo_block,"min","max", "6B call_coarse");

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  Solver * solve_coarse = simulation->problem()->solver(index_solve_coarse_);

  solve_coarse->set_min_level(min_level_coarse_);
  solve_coarse->set_max_level(max_level_coarse_);
  solve_coarse->set_sync_id (enzo_sync_id_solver_mg0_coarse);
  solve_coarse->set_callback(CkIndex_EnzoBlock::p_solver_mg0_solve_coarse());
  
  //  TRACE_FIELD_("X",X,1.0);
  //  TRACE_FIELD_("B",B,1.0);
  DEBUG_FIELD (enzo_block,ib_,"B solve_coarse");
  solve_coarse->apply(A_,ix_,ib_,enzo_block);
  DEBUG_FIELD (enzo_block,ix_,"X solve_coarse");
}

//----------------------------------------------------------------------

void EnzoSolverMg0::restrict_send(EnzoBlock * enzo_block) throw()
///
///      A.residual(R,B,X)
///      pack R
///      index_parent.p_restrict_recv(R)
{
  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "7 restrict_send");

  TRACE_LEVEL("EnzoSolverMg0::restrict_send",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_send()");

  Data * data = enzo_block->data();
  Field field = data->field();

  A_->residual(ir_, ib_, ix_, enzo_block);
  DEBUG_FIELD (enzo_block,ib_,"B restrict_send");
  DEBUG_FIELD (enzo_block,ir_,"R restrict_send");
  DEBUG_FIELD (enzo_block,ix_,"X restrict_send");

  if ( is_finest_(enzo_block) ) {
    enzo_float * R = (enzo_float*) field.values(ir_);
    for (int iz=gz_; iz<nz_+gz_; iz++) {
      for (int iy=gy_; iy<ny_+gy_; iy++) {
	for (int ix=gx_; ix<nx_+gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  rr_local_ += R[i]*R[i];
	}
      }
    }
  }

#ifdef DEBUG_COPY
  enzo_float * R_copy = (enzo_float*) field.values("RMG");
  enzo_float * R = (enzo_float*) field.values(ir_);
  double rsum=0.0;
  for (int i=0; i<mx_*my_*mz_; i++) {
    R_copy[i]=R[i];
    rsum += std::abs(R[i]);
  }
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("DEBUG_COPY rsum = %g\n",rsum);
#endif

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
#ifdef DEBUG_FIELD_FACE  
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%d %s DEBUG_FIELD_FACE creating %p\n",CkMyPe(),__FILE__,field_face);
#endif

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
{
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY enter p_solver_mg0_restrict_recv\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  DEBUG_FIELD_MSG(this,"restrict_inter");

  TRACE_MG_BLOCK(this,"EnzoBlock::restrict_recv()");
  DEBUG_INDEX(this,"restrict_inter",NULL);

  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  solver->restrict_recv(this,msg);

#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  exit p_solver_mg0_restrict_recv\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  performance_stop_(perf_compute,__FILE__,__LINE__);

}

//----------------------------------------------------------------------

void EnzoSolverMg0::restrict_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
/// 
///      [ unpack B ]
///      if (sync.next())
///          begin_cycle()
{

  SOLVER_CONTROL(enzo_block,"coarse","fine-1", "8 restrict_recv");

  DEBUG_FIELD_MSG(enzo_block,"restrict_recv");

  // Unpack "B" vector data from children

  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  Refresh * refresh = new Refresh;
  refresh->add_field(ib_);

  DEBUG_FIELD (enzo_block,ib_,"B restrict_recv");

  // copy data from msg to this EnzoBlock

  int * ic3 = msg->ic3;

  FieldFace * field_face = enzo_block->create_face 
    (if3, ic3, lg3, refresh_coarse, refresh, true);
#ifdef DEBUG_FIELD_FACE  
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%d %s DEBUG_FIELD_FACE creating %p\n",CkMyPe(),__FILE__,field_face);
#endif

  field_face->set_restrict(restrict());

  char * a = msg->a;
  field_face->array_to_face(a, enzo_block->data()->field());
  delete field_face;

  delete msg;

  // continue with EnzoSolverMg0

  TRACE_LEVEL("EnzoSolverMg0::restrict_recv",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_recv()");

  //  TRACE_FIELD_("B",B,1.0);

#ifdef NEW_SYNC  
  ScalarData<Sync> * scalar_data  = enzo_block->data()->scalar_data_sync();
  ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
  Sync & sync = scalar_data->value(scalar_descr,isync_restrict_);
  if (sync.next())
    {
      if (AFTER_CYCLE(enzo_block,CYCLE)) CkPrintf ("DEBUG_BEGIN_CYCLE %s %s 2\n",name_.c_str(),enzo_block->name().c_str());
      begin_cycle_ (enzo_block);
    }
#else
  if (enzo_block->mg_sync_restrict_next()) 
    {
      if (AFTER_CYCLE(enzo_block,CYCLE)) CkPrintf ("DEBUG_BEGIN_CYCLE %s %s 3\n",name_.c_str(),enzo_block->name().c_str());
      begin_cycle_ (enzo_block);
    }
#endif

}

//----------------------------------------------------------------------

void EnzoSolverMg0::solve_coarse(EnzoBlock * enzo_block) throw()
/// 
///      solve A X = B
///      end_cycle()
{

  SOLVER_CONTROL(enzo_block,"min","max", "9 solve_coarse");

#ifdef DEBUG_COPY
  Field field = enzo_block->data()->field();
  enzo_float * B_copy = (enzo_float*) field.values("BMG");
  enzo_float * B = (enzo_float*) field.values(ib_);
  double bsum=0.0;
  for (int i=0; i<mx_*my_*mz_; i++) {
    B_copy[i]=B[i];
    bsum += std::abs(B[i]);
  }
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("DEBUG_COPY bsum = %g\n",bsum);
#endif
  
#ifdef DEBUG_COPY
  enzo_float * X_copy = (enzo_float*) field.values("XMG");
  enzo_float * X = (enzo_float*) field.values(ix_);
  double xsum=0.0;
  for (int i=0; i<mx_*my_*mz_; i++) {
    X_copy[i]=X[i];
    xsum += std::abs(X[i]);
  }
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("DEBUG_COPY xsum = %g\n",xsum);
#endif
  
  TRACE_LEVEL("EnzoSolverMg0::solve_coarse",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::solver_coarse()");
 
  /// Prolong solution to next-finer level

  const int level = enzo_block->level();

  if (level == min_level_) {

    if ( ! is_finest_(enzo_block) ) {

#ifdef DEBUG_PROLONG
      if (AFTER_CYCLE(enzo_block,CYCLE))
	CkPrintf ("%s %s %s DEBUG_PROLONG call prolong_send_()\n",
		  __FILE__,enzo_block->name().c_str(),name_.c_str());
#endif  
      prolong_send_ (enzo_block);
      
    }

    end_cycle (enzo_block);

  } else if (level > min_level_ && is_active_(enzo_block)) {

#ifdef DEBUG_PROLONG
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%s %s %s DEBUG_PROLONG A call prolong_recv(NULL)\n",
		__FILE__,enzo_block->name().c_str(),name_.c_str());
#endif  
    enzo_block->solver_mg0_prolong_recv(NULL);
  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::prolong_send_(EnzoBlock * enzo_block) throw()
/// 
///      for child           
///         pack X
///         child.prolong_recv(X)
{
#ifdef DEBUG_PROLONG
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%s %s %s DEBUG_PROLONG enter prolong_send()\n",
	      __FILE__,enzo_block->name().c_str(),name_.c_str());
#endif  
  TRACE_LEVEL("EnzoSolverMg0::prolong_send",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::prolong_send()");

  SOLVER_CONTROL(enzo_block,"coarse","fine-1", "10 prolong_send");

  ItChild it_child(enzo_block->rank());
  int ic3[3];
  //  TRACE_FIELD_("X",X,1.0);
  DEBUG_FIELD (enzo_block,ix_,"X prolong_send");
  while (it_child.next(ic3)) {

    Index index_child = enzo_block->index().index_child(ic3,min_level_);

    // Pack and send "X" to children

    // <COMMON CODE> in restrict_send_() and prolong_send_()

    int if3[3] = {0,0,0};
    bool lg3[3] = {true,true,true};
    Refresh * refresh = new Refresh;
    refresh->add_field(ix_);
    
    // copy data from EnzoBlock to array via FieldFace

    FieldFace * field_face = enzo_block->create_face
      (if3, ic3, lg3, refresh_fine, refresh, true);
#ifdef DEBUG_FIELD_FACE  
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%d %s DEBUG_FIELD_FACE creating %p\n",CkMyPe(),__FILE__,field_face);
#endif

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

#ifdef DEBUG_PROLONG
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%s %s %s DEBUG_PROLONG B call prolong_recv()\n",
		__FILE__,enzo_block->name().c_str(),name_.c_str());
#endif  
    enzo_block->thisProxy[index_child].p_solver_mg0_prolong_recv(msg);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  solver_mg0_prolong_recv(msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
  
}

//----------------------------------------------------------------------

void EnzoBlock::solver_mg0_prolong_recv(FieldMsg * msg)
{
#ifdef DEBUG_PROLONG
#ifdef NEW_SYNC
#else
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%s %s %s DEBUG_PROLONG enter p_solver_mg0_prolong_recv(%d/%d)\n",
	      __FILE__,name().c_str(),name_.c_str(),
	      mg_sync_prolong_.value(),mg_sync_prolong_.stop(),name_.c_str());
#endif
#endif  
  TRACE_MG (this,"EnzoBlock::p_solver_mg0_prolong_recv()");


#ifdef NEW_SYNC
#else  
  // Save message
  if (msg != NULL) mg_msg_ = msg;

  // Return if not ready yet
  if (! mg_sync_prolong_next()) return;

  // Restore saved message
  msg = mg_msg_;
  mg_msg_ = NULL;
#endif  
  
  static_cast<EnzoSolverMg0*> (solver())->prolong_recv(this,msg);

#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  exit p_solver_mg0_prolong_recv\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
}

//----------------------------------------------------------------------

void EnzoSolverMg0::prolong_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
/// 
///      [ unpack C ]
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{

  // Save message

  // Return if not ready yet
#ifdef NEW_SYNC
  if (msg != NULL) enzo_block->mg_set_msg (msg);
  ScalarData<Sync> * scalar_data = enzo_block->data()->scalar_data_sync();
  ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
  Sync & sync = scalar_data->value(scalar_descr,isync_prolong_);
  if (! sync.next() ) return;
  // Restore saved message
  msg = enzo_block->mg_get_msg();
  enzo_block->mg_set_msg(NULL);
#endif

  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "11 prolong_recv");
  
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%d %s %p %s DEBUG_ENTRY enter p_solver_mg0_prolong_recv\n",
	      CkMyPe(),name().c_str(),this,name_.c_str());
#endif

  // Unpack "C" vector data from children

#ifdef DEBUG_PROLONG
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%s %s %s DEBUG_PROLONG enter prolong_recv()\n",
	      __FILE__,enzo_block->name().c_str(),name_.c_str());
#endif  
  DEBUG_FIELD_MSG(enzo_block,"prolong_recv");
    
  int if3[3] = {0,0,0};
  bool lg3[3] = {true,true,true};
  Refresh * refresh = new Refresh;
  refresh->add_field(ic_);

  // copy data from msg to this EnzoBlock

  FieldFace * field_face = enzo_block->create_face 
    (if3, msg->ic3, lg3, refresh_fine, refresh, true);
#ifdef DEBUG_FIELD_FACE  
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%d %s DEBUG_FIELD_FACE creating %p\n",CkMyPe(),__FILE__,field_face);
#endif

  field_face->set_prolong(prolong());

  field_face->array_to_face (msg->a, enzo_block->data()->field());

  delete field_face;

  // GET DATA
  
  delete msg;

  TRACE_LEVEL("EnzoSolverMg0::prolong_recv",enzo_block);
  TRACE_MG (enzo_block,"EnzoSolverMg0::prolong_recv()");
  
  Field field = enzo_block->data()->field();

  enzo_float * X = (enzo_float*) field.values(ix_);
  enzo_float * C = (enzo_float*) field.values(ic_);

#ifdef DEBUG_COPY
  enzo_float * C_copy = (enzo_float*) field.values("CMG");
  double csum=0.0;
  for (int i=0; i<mx_*my_*mz_; i++) {
    C_copy[i]=C[i];
    csum += std::abs(C[i]);
  }
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("DEBUG_COPY csum = %g\n",csum);
 
#endif
  TRACE_FIELD_("C",C,1.0);

  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] += C[i];
      }
    }
  }

  DEBUG_FIELD (enzo_block,ic_,"C update");
  DEBUG_FIELD (enzo_block,ix_,"X update");
  if (index_smooth_post_ >= 0) {
    Simulation * simulation = proxy_simulation.ckLocalBranch();
    Solver * smooth_post = simulation->problem()->solver(index_smooth_post_);

    smooth_post->set_min_level(enzo_block->level());
    smooth_post->set_max_level(enzo_block->level());
    smooth_post->set_sync_id (enzo_sync_id_solver_mg0_post);
    smooth_post->set_callback(CkIndex_EnzoBlock::p_solver_mg0_post_smooth());
  
    smooth_post->apply(A_,ix_,ib_,enzo_block);

  } else {

    post_smooth (enzo_block);
  }

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_post_smooth()
{
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY enter p_solver_mg0_post_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  TRACE_MG (this,"EnzoBlock::p_solver_mg0_post_smooth()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  solver->post_smooth(enzo_block);
  
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  exit p_solver_mg0_post_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::post_smooth(EnzoBlock * enzo_block) throw()
///
///      smooth.apply (A,X,B)
///      end_cycle()
{
  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "12 post-smooth");

  TRACE_LEVEL("EnzoSolverMg0::post_smooth",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::post_smooth()");

  const int level = enzo_block->level();

  if ( ! is_finest_(enzo_block) ) {

    prolong_send_ (enzo_block);
  } 

  TRACE_MG(enzo_block,"EnzoSolverMg0::post_smooth() calling end_cycle()");
  end_cycle (enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::end_cycle(EnzoBlock * enzo_block) throw()
/// 
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
{
  SOLVER_CONTROL(enzo_block,"min","max", "13 end_cycle");
  TRACE_LEVEL("EnzoSolverMg0::end_cycle",enzo_block);

  TRACE_MG(enzo_block,"EnzoSolverMg0::end_cycle()");
  
  enzo_block->mg_iter_increment();

  bool is_converged = is_converged_(enzo_block);
  bool is_diverged  = is_diverged_(enzo_block);

  const int iter = enzo_block->mg_iter();
	    
  const int level = enzo_block->level();

  const bool l_output =
    ( ( enzo_block->index().is_zero() && level == max_level_) &&
      ( (is_converged) || (is_diverged) ||
	(monitor_iter_ && (iter % monitor_iter_) == 0 )) );

  if (l_output) {
    monitor_output_(enzo_block,iter,rr0_,0.0,rr_,0.0);
  }

  if (is_converged || is_diverged) {

    // Do an optional final smoothing on the full mesh For use in Dan
    // Reynolds HG algorithm in which Mg0 with no pre- or
    // post-smoothings is used as a preconditioner to BiCgStab
    
    if (index_smooth_last_ >= 0 && (is_finest_(enzo_block)) ) {

      Simulation * simulation = proxy_simulation.ckLocalBranch();
      Solver * smooth_last = simulation->problem()->solver(index_smooth_last_);
      smooth_last->set_sync_id (enzo_sync_id_solver_mg0_last);
      smooth_last->set_callback(CkIndex_EnzoBlock::p_solver_mg0_last_smooth());
  
      smooth_last->apply(A_,ix_,ib_,enzo_block);

    } else {

      end (enzo_block);
    }

  } else if ( is_finest_(enzo_block)) {

    begin_cycle_ (enzo_block);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_last_smooth()
{
  performance_start_(perf_compute,__FILE__,__LINE__);
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY enter p_solver_mg0_last_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  TRACE_MG (this,"EnzoBlock::p_solver_mg0_last_smooth()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);

  solver->end(enzo_block);
  
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  exit p_solver_mg0_last_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  performance_stop_(perf_compute,__FILE__,__LINE__);
}


//======================================================================

bool EnzoSolverMg0::is_converged_(EnzoBlock * enzo_block) const
{
  TRACE_MG(enzo_block,"EnzoSolverMg0::is_converged");
  return (rr0_ != 0.0 && rr_/rr0_ < res_tol_);
}

//----------------------------------------------------------------------

bool EnzoSolverMg0::is_diverged_(EnzoBlock * enzo_block) const
/// [*]
{
  TRACE_MG(enzo_block,"EnzoSolverMg0::is_diverged");
  return (enzo_block->mg_iter() >= iter_max_);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::end(Block * block)
{
  SOLVER_CONTROL(block,"min","max", "14 end");

  TRACE_MG(block,"EnzoSolverMg0::end");
    
  Field field = block->data()->field();
  
  deallocate_temporary_(field,block);
    
  Solver::end_(block);

  CkCallback(callback_,
	     CkArrayIndexIndex(block->index()),
	     block->proxy_array()).send();

}
