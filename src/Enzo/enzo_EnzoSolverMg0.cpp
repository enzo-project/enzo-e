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

#include "cello.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

// #define DEBUG_BARRIER
// #define DEBUG_SCALARS
// #define DEBUG_COPY

// #define DEBUG_ENTRY
// #define DEBUG_SOLVER_MG0
// #define DEBUG_SOLVER_CONTROL
// #define DEBUG_TRACE_LEVEL
// #define DEBUG_FIELD_MESSAGE
// #define DEBUG_PROLONG


#define AFTER_CYCLE(BLOCK,CYCLE) (BLOCK->cycle() >= CYCLE)
#define CYCLE 000

#ifdef DEBUG_SCALARS
#   define TRACE_SCALARS(MSG,BLOCK)					\
  {									\
    if (AFTER_CYCLE(BLOCK,CYCLE))					\
      CkPrintf ("%s %s DEBUG_SCALARS restrict %d/%d prolong %d/%d iter %d rr %lg rr_local %lg rr0 %lg %s\n", \
		name_.c_str(),MSG,					\
		psync_restrict(BLOCK)->value(),				\
		psync_restrict(BLOCK)->stop(),				\
		psync_prolong(BLOCK)->value(),				\
		psync_prolong(BLOCK)->stop(),				\
		*piter(BLOCK), rr_,rr_local_,rr0_,			\
		BLOCK->name().c_str()    );				\
  }
# else
#   define TRACE_SCALARS(MSG,BLOCK) /* ... */
#endif

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

#define DEBUG_FIELD(BLOCK,ID,NAME)			\
  {							\
    if (AFTER_CYCLE(BLOCK,CYCLE)) {			\
      Data * data = BLOCK->data();			\
      Field field = data->field();			\
      enzo_float * X = (enzo_float*) field.values(ID);	\
      double sum_a=0.0;						\
      double sum_abs=0.0;						\
      for (int iz=gz_; iz<mz_-gz_; iz++) {				\
	for (int iy=gy_; iy<my_-gy_; iy++) {				\
	  for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	    int i = ix + mx_*(iy + my_*iz);				\
	    sum_a+=X[i];						\
	    sum_abs+=std::abs(X[i]);					\
	  }								\
	}								\
      }									\
      CkPrintf ("%s:%d %s %s COPY_FIELD %d %s shift %20.15lg %20.15lg\n" \
		,__FILE__,__LINE__,BLOCK->name().c_str(),name().c_str(),ID,NAME,sum_a, sum_abs); \
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


#ifdef DEBUG_SOLVER_CONTROL
#   define SOLVER_CONTROL(BLOCK,MIN,MAX,MESSAGE)			\
  if (AFTER_CYCLE(BLOCK,CYCLE))						\
    CkPrintf ("DEBUG_SOLVER_CONTROL %-4s %-10s %s\n",	\
	      name_.c_str(),BLOCK->name().c_str(),MESSAGE);
#else
#   define SOLVER_CONTROL(BLOCK,MIN,MAX,MESSAGE) /* ... */
#endif

#ifdef DEBUG_COPY
#  define DEBUG_COPY_FIELD(NAME,ID)					\
  if (AFTER_CYCLE(enzo_block,CYCLE))					\
    {									\
      Field field = enzo_block->data()->field();			\
      char name[20];							\
      sprintf(name,"%s%s%d",this->name().c_str(), NAME,enzo_block->level()-min_level_); \
      CkPrintf ("%s Copying field %d to %s\n",enzo_block->name().c_str(),ID,name); \
      enzo_float * VALUES_copy = (enzo_float*) field.values(name);	\
      enzo_float * VALUES = (enzo_float*) field.values(ID);		\
      double sum_a=0.0;						\
      double sum_abs=0.0;						\
      for (int iz=gz_; iz<mz_-gz_; iz++) {				\
      for (int iy=gy_; iy<my_-gy_; iy++) {				\
	for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	  int i=ix+mx_*(iy+my_*iz);					\
	  VALUES_copy[i]=VALUES[i];					\
	  sum_a+=X[i];							\
	  sum_abs+=std::abs(X[i]);					\
      }									\
	CkPrintf ("%s:%d %s COPY_FIELD %d %s shift %20.15lg\n"		\
		  ,__FILE__,__LINE__,BLOCK->name().c_str(),ID,COPY, \
		  sum_abs!=0?sum_a/sum_abs:0.0);		    \
    }
#else
#  define DEBUG_COPY_FIELD(NAME,ID) /* ... */
#endif

#ifdef DEBUG_BARRIER
#   define TRACE_BARRIER(BLOCK,SOLVER,MSG) \
  CkPrintf ("DEBUG_BARRIER %s %s %s\n", \
	    BLOCK->name().c_str(),SOLVER->name().c_str(),MSG);
#else
#   define TRACE_BARRIER(BLOCK,SOLVER,MSG)  /* ... */
#endif


//======================================================================

EnzoSolverMg0::EnzoSolverMg0
(std::string name,
 std::string field_x,
 std::string field_b,
 int monitor_iter,
 int restart_cycle,
 int solve_type,
 int min_level,
 int max_level,
 int iter_max,
 double res_tol,
 int index_smooth_pre,
 int index_solve_coarse,
 int index_smooth_post,
 int index_smooth_last,
 Restrict * restrict,
 Prolong * prolong,
 int coarse_level) 
  : Solver(name,
	   field_x,
	   field_b,
	   monitor_iter,
	   restart_cycle,
	   solve_type,
	   min_level,
	   max_level),
    bs_(0), bc_(0),
    rr_(0), rr_local_(0), rr0_(0),
    res_tol_(res_tol),
    A_(NULL),
    index_smooth_pre_(index_smooth_pre),
    index_solve_coarse_(index_solve_coarse),
    index_smooth_post_(index_smooth_post),
    index_smooth_last_(index_smooth_last),
    restrict_(restrict),
    prolong_(prolong),
    iter_max_(iter_max), 
    ic_(-1), ir_(-1),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    coarse_level_(coarse_level)
{
  // Initialize temporary fields

  FieldDescr * field_descr = cello::field_descr();

  ir_ = field_descr->insert_temporary();
  ic_ = field_descr->insert_temporary();

  /// Initialize default Refresh

  const int ir = add_refresh
    (4,0,neighbor_type_(),sync_barrier,
     enzo_sync_id_solver_mg0);

  refresh(ir)->add_field (ix_);  // NOTE: ix_ set in Solver::Solver()
  refresh(ir)->add_field (ir_);
  refresh(ir)->add_field (ic_);

  ScalarDescr * scalar_descr_int  = cello::scalar_descr_int();
  i_iter_  = scalar_descr_int ->new_value(name + ":iter");
  
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  i_sync_restrict_ = scalar_descr_sync->new_value(name + ":restrict");
  i_sync_prolong_  = scalar_descr_sync->new_value(name + ":prolong");

  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  i_msg_ = scalar_descr_void->new_value(name + ":msg");

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

void EnzoSolverMg0::apply ( std::shared_ptr<Matrix> A, Block * block) throw()
{
  TRACE_MG(block,"EnzoSolverMg0::apply()");

  Solver::begin_(block);

  A_ = A;

  allocate_temporary_(block);

  // clear scalars
  bs_ = 0.0;
  bc_ = 0.0;
  rr_ = 0.0;
  rr_local_ = 0.0;
  rr0_ = 0.0;
  *piter(block) = 0.0;
  *pmsg(block) = NULL;

  /// Current and initial residual norm R'*R

  TRACE_LEVEL("EnzoSolverMg0::apply",block);

  Field field = block->data()->field();

  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  EnzoBlock * enzo_block = enzo::block(block);

  // Initialize sync counters for restrict and prolong
  
  Sync * sync_restrict = psync_restrict(block);

  sync_restrict->reset();
  sync_restrict->set_stop(cello::num_children());
  
  Sync * sync_prolong = psync_prolong(block);

  sync_prolong->reset();
  sync_prolong->set_stop(2); // self and parent

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

  *piter(enzo_block) = 0.0;

  Field field = enzo_block->data()->field();

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

      compute_shift_(enzo_block,reduce);
    }

#ifdef DEBUG_SOLVER_MG0    
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%s DEBUG_SOLVER_MG0 %s bs %Lf bc %Lf\n",
		enzo_block->name().c_str(),name_.c_str(),reduce[0],reduce[1]);
#endif    

    /// initiate callback for p_solver_begin_solve and contribute to
    /// sum and count

    CkCallback callback(CkIndex_EnzoBlock::r_solver_mg0_begin_solve(NULL), 
			enzo::block_array());

    SOLVER_CONTROL (enzo_block,"min","max","1 calling begin_solve (shift)");

    TRACE_BARRIER(enzo_block,this,"shift");
    
    enzo_block->contribute(2*sizeof(long double), &reduce, 
			   sum_long_double_2_type, callback);

  } else {

    SOLVER_CONTROL(enzo_block,"min","max","2 calling begin_solve (no shift)");

    begin_solve (enzo_block,NULL);

  }

}

//----------------------------------------------------------------------

void EnzoSolverMg0::compute_shift_
(EnzoBlock * enzo_block,long double * reduce) throw()
{
  Field field = enzo_block->data()->field();

  enzo_float* B = (enzo_float*) field.values(ib_);

  for (int iz=gz_; iz<mz_-gz_; iz++) {
    for (int iy=gy_; iy<my_-gy_; iy++) {
      for (int ix=gx_; ix<mx_-gx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	reduce[0] += B[i];
	reduce[1] += 1.0;
      }
    }
  }
}
//----------------------------------------------------------------------

void EnzoBlock::r_solver_mg0_begin_solve(CkReductionMsg* msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);

  static_cast<EnzoSolverMg0*> (solver())->begin_solve(this,msg);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::begin_solve(EnzoBlock * enzo_block,
				CkReductionMsg *msg) throw()
{
  DEBUG_FIELD(enzo_block,ix_,"X begin_solve");
  DEBUG_FIELD(enzo_block,ib_,"B begin_solve");

  TRACE_SCALARS("begin_solve",enzo_block);
  TRACE_LEVEL("EnzoSolverMg0::begin_solve",enzo_block);
  // start the MG V-cycle with the max level blocks


  SOLVER_CONTROL(enzo_block,"min","max", "3 calling do shift");
  
  do_shift_(enzo_block,msg);

  // control flow starts at leaves, even in level > max_level,
  // since coarse solve may require reductions over all Blocks

  if (is_finest_(enzo_block)) {

    SOLVER_CONTROL(enzo_block, "fine","max", "4 calling begin_cycle");

    begin_cycle_ (enzo_block);

  } else {

    int level = enzo_block->level();
    if ( ! (coarse_level_ <= level && level <= max_level_) ) {
      DEBUG_FIELD (enzo_block,ix_,"X begin_solve");
      call_coarse_solver(enzo_block);
    }

  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::do_shift_(EnzoBlock * enzo_block,
			      CkReductionMsg *msg) throw()
{
  if (msg != NULL) {
    
    long double* data = (long double*) msg->getData();

    bs_ = data[0];
    bc_ = data[1];

    delete msg;
  } 
  
  if (A_->is_singular() && is_finest_(enzo_block)) {

    // Shift B if needed to be in range(A) for periodic b.c.
    
    SOLVER_CONTROL(enzo_block,"fine","fine", "5 applying shift");

    Field field = enzo_block->data()->field();

    double shift = -bs_ / bc_;
#ifdef DEBUG_SOLVER_MG0      
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%s DEBUG_SOLVER_MG0 bs B %lf bc %lf\n",
		enzo_block->name().c_str(),double(bs_),double(bc_));
#endif      

    DEBUG_FIELD (enzo_block,ib_,"B shift");
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
  TRACE_SCALARS("begin_cycle",enzo_block);
  monitor_output_(enzo_block);
  
  TRACE_LEVEL("EnzoSolverMg0::begin_cycle",enzo_block);

  TRACE_MG(enzo_block,"EnzoSolverMg0::begin_cycle()");

  Field field = enzo_block->data()->field();
  
  const int level = enzo_block->level();

  if ( ! is_finest_(enzo_block) ) {
    enzo_float * X = (enzo_float*) field.values(ix_);
    std::fill_n(X,mx_*my_*mz_,0.0);
  }

  if (level == coarse_level_) {

    SOLVER_CONTROL(enzo_block,"min","coarse","6 calling coarse solve");

    TRACE_MG(enzo_block,"calling coarse solve");

    DEBUG_FIELD(enzo_block,ix_,"X clear");

    call_coarse_solver(enzo_block);

  } else {

    TRACE_MG(enzo_block,"calling smoother");

    // if ( ! is_finest_(enzo_block) ) {

    //   enzo_float * X = (enzo_float*) field.values(ix_);
    //   std::fill_n(X,mx_*my_*mz_,0.0);
    //   DEBUG_FIELD(enzo_block,ix_,"X clear");

    // }

#ifdef DEBUG_SOLVER_REFRESH    
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("DEBUG_SOLVER_REFRESH refresh sync_face %d\n",refresh.sync_id());
#endif

    if (index_smooth_pre_ >= 0) {

      SOLVER_CONTROL(enzo_block,"coarse+1","fine", "7 calling pre-smooth");
      call_pre_smoother (enzo_block);

    } else {

      SOLVER_CONTROL(enzo_block,"coarse+1","fine", "8 calling restrict");
      restrict (enzo_block);

    }

  }
}

//----------------------------------------------------------------------

void EnzoSolverMg0::monitor_output_(EnzoBlock * enzo_block)
{
  const int iter = *(piter(enzo_block));

  const int level = enzo_block->level();

  const bool l_output =
    ( ( enzo_block->index().is_root()) &&
      ( (iter == 0))); // ||

  if (l_output) {
    Solver::monitor_output_(enzo_block,iter,rr0_,0.0,rr_,0.0);
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
    CkPrintf ("DEBUG_SOLVER_MG0 %s %s solver->rr_local() = %lg\n",
	      name().c_str(),solver->name().c_str(),solver->rr_local());
#endif

  CkCallback callback(CkIndex_EnzoBlock::r_solver_mg0_barrier(NULL), 
		      enzo::block_array());
  long double data[1] = {solver->rr_local()};
  TRACE_BARRIER(this,solver,"barrier");
  contribute(sizeof(long double), data,  sum_long_double_type, callback);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_mg0_barrier(CkReductionMsg* msg)
{
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

#ifdef DEBUG_SOLVER_MG0
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("DEBUG_SOLVER_MG0 %s %s rr_ %lg\n",
	      name().c_str(),solver->name().c_str(),solver->rr());
#endif

  performance_start_(perf_compute,__FILE__,__LINE__);

  long double rr = ((long double*) msg->getData())[0];
  solver->set_rr(rr);
  solver->set_rr_local(0.0);
  if (*solver->piter(this)==0) solver->set_rr0(rr);
  
  delete msg;

#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  after barrier\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  TRACE_MG_BLOCK(this,"EnzoBlock::solver_mg0_coarse_solve()");
  
  solver->prolong(this);

  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_restrict()
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

  solver->restrict(this);

#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  exit p_solver_mg0_pre_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::restrict(EnzoBlock * enzo_block) throw()
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "9 restrict");
  
  TRACE_LEVEL("EnzoSolverMg0::restrict",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict()");

  restrict_send (enzo_block);

  // All Blocks must call coarse solver since may involve
  // global reductions

  DEBUG_FIELD (enzo_block,ix_,"X restrict");
  call_coarse_solver(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_coarse_solver(EnzoBlock * enzo_block) throw()
{
  SOLVER_CONTROL(enzo_block,"min","max", "10 call_coarse_solver");

  Solver * solve_coarse = cello::solver(index_solve_coarse_);

  solve_coarse->set_min_level(min_level_);
  solve_coarse->set_max_level(coarse_level_);
  solve_coarse->set_sync_id (enzo_sync_id_solver_mg0_coarse);
  solve_coarse->set_callback(CkIndex_EnzoBlock::p_solver_mg0_solve_coarse());
  
  DEBUG_FIELD (enzo_block,ib_,"B solve_coarse");
  DEBUG_FIELD (enzo_block,ix_,"X solve_coarse 1");

  solve_coarse->set_field_x (ix_);
  solve_coarse->set_field_b (ib_);

  solve_coarse->apply(A_,enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_pre_smoother(EnzoBlock * enzo_block) throw()
{
  SOLVER_CONTROL(enzo_block,"min","max", "11 call_pre_smoother");
  Solver * smooth_pre = cello::solver(index_smooth_pre_);

  smooth_pre->set_min_level(enzo_block->level());
  smooth_pre->set_max_level(enzo_block->level());
  smooth_pre->set_sync_id (enzo_sync_id_solver_mg0_pre);
  smooth_pre->set_callback(CkIndex_EnzoBlock::p_solver_mg0_restrict());

  DEBUG_FIELD (enzo_block,ib_,"B pre_smooth");

  smooth_pre->set_field_x(ix_);
  smooth_pre->set_field_b(ib_);
  
  smooth_pre->apply(A_,enzo_block);
  
  DEBUG_FIELD (enzo_block,ix_,"X pre_smooth");
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_post_smoother(EnzoBlock * enzo_block) throw()
{
  SOLVER_CONTROL(enzo_block,"min","max", "12 call_post_smooth");
  Solver * smooth_post = cello::solver(index_smooth_post_);

  smooth_post->set_min_level(enzo_block->level());
  smooth_post->set_max_level(enzo_block->level());
  smooth_post->set_sync_id (enzo_sync_id_solver_mg0_post);
  smooth_post->set_callback(CkIndex_EnzoBlock::p_solver_mg0_post_smooth());

  DEBUG_FIELD (enzo_block,ib_,"B post_smooth");

  smooth_post->set_field_x(ix_);
  smooth_post->set_field_b(ib_);
  
  smooth_post->apply(A_,enzo_block);
  
  DEBUG_FIELD (enzo_block,ix_,"X post_smooth");
}

//----------------------------------------------------------------------

void EnzoSolverMg0::call_last_smoother(EnzoBlock * enzo_block) throw()
{
  SOLVER_CONTROL(enzo_block,"min","max", "13 call_last_smooth");
  Solver * smooth_last = cello::solver(index_smooth_last_);

  smooth_last->set_sync_id (enzo_sync_id_solver_mg0_last);
  smooth_last->set_callback(CkIndex_EnzoBlock::p_solver_mg0_last_smooth());

  DEBUG_FIELD (enzo_block,ib_,"B last_smooth");

  smooth_last->set_field_x(ix_);
  smooth_last->set_field_b(ib_);

  smooth_last->apply(A_,enzo_block);

  DEBUG_FIELD (enzo_block,ix_,"X last_smooth");
}

//----------------------------------------------------------------------

void EnzoSolverMg0::restrict_send(EnzoBlock * enzo_block) throw()
///
///      A.residual(R,B,X)
///      pack R
///      index_parent.p_restrict_recv(R)
{

  compute_residual_(enzo_block);

  DEBUG_FIELD(enzo_block,ir_,"R restrict_send");
  FieldMsg * msg = pack_residual_(enzo_block);
  
  DEBUG_FIELD_MSG(enzo_block,"restrict_send");

  Index index_parent = enzo_block->index().index_parent(min_level_);

  enzo::block_array()[index_parent].p_solver_mg0_restrict_recv(msg);

}

//----------------------------------------------------------------------

void EnzoSolverMg0::compute_residual_(EnzoBlock * enzo_block) throw()
{
  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "14 compute_residual");

  TRACE_LEVEL("EnzoSolverMg0::restrict_send",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_send()");

  Field field = enzo_block->data()->field();

  A_->residual(ir_, ib_, ix_, enzo_block);

  DEBUG_FIELD (enzo_block,ir_,"R residual");

  if ( is_finest_(enzo_block) ) {
    enzo_float * R = (enzo_float*) field.values(ir_);
    for (int iz=gz_; iz<mz_-gz_; iz++) {
      for (int iy=gy_; iy<my_-gy_; iy++) {
	for (int ix=gx_; ix<mx_-gx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  rr_local_ += R[i]*R[i];
	}
      }
    }
  }

  DEBUG_COPY_FIELD("R",ir_);
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

  SOLVER_CONTROL(enzo_block,"coarse","fine-1", "15 restrict_recv");

  DEBUG_FIELD_MSG(enzo_block,"restrict_recv");

  // Unpack "B" vector data from children

  DEBUG_FIELD (enzo_block,ix_,"X update");
  unpack_residual_(enzo_block,msg);
  DEBUG_FIELD (enzo_block,ix_,"X update");
  
  DEBUG_FIELD(enzo_block,ir_,"R restrict_recv");

  // continue with EnzoSolverMg0

  TRACE_LEVEL("EnzoSolverMg0::restrict_recv",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_recv()");

  if (psync_restrict(enzo_block)->next())
    {
      SOLVER_CONTROL(enzo_block,"coarse","fine-1", "16 call begin_cycle");
      begin_cycle_ (enzo_block);
    }

}

//----------------------------------------------------------------------

void EnzoSolverMg0::prolong(EnzoBlock * enzo_block) throw()
/// 
///      solve A X = B
///      end_cycle()
{

  SOLVER_CONTROL(enzo_block,"min","max", "17 prolong");
  DEBUG_FIELD (enzo_block,ix_,"X solve_coarse 2");

  DEBUG_COPY_FIELD("B",ib_);
  DEBUG_COPY_FIELD("X",ix_);
  
  TRACE_LEVEL("EnzoSolverMg0::solve_coarse",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::solver_coarse()");
 
  /// Prolong solution to next-finer level

  const int level = enzo_block->level();

  if (level == coarse_level_) {

    if ( ! is_finest_(enzo_block) ) {

      SOLVER_CONTROL(enzo_block,"coarse","fine-1", "18 call prolong_send");
#ifdef DEBUG_PROLONG
      if (AFTER_CYCLE(enzo_block,CYCLE))
	CkPrintf ("%s %s %s DEBUG_PROLONG call prolong_send_()\n",
		  __FILE__,enzo_block->name().c_str(),name_.c_str());
#endif  
      prolong_send_ (enzo_block);
      
    }
  }

  if (coarse_level_ < level && level <= max_level_) {
#ifdef DEBUG_PROLONG
    if (AFTER_CYCLE(enzo_block,CYCLE))
      CkPrintf ("%s %s %s DEBUG_PROLONG A call prolong_recv(NULL)\n",
		__FILE__,enzo_block->name().c_str(),name_.c_str());
#endif  
    SOLVER_CONTROL(enzo_block,"coarse","fine-1", "20 call prolong_recv");
    enzo_block->solver_mg0_prolong_recv(NULL);
  } else {
  
    SOLVER_CONTROL(enzo_block,"coarse","fine-1", "19 call end_cycle");
    end_cycle (enzo_block);
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

  SOLVER_CONTROL(enzo_block,"coarse","fine-1", "22 prolong_send");

  ItChild it_child(cello::rank());
  int ic3[3];
  DEBUG_FIELD (enzo_block,ix_,"X prolong_send");
  while (it_child.next(ic3)) {

    FieldMsg * msg = pack_correction_(enzo_block,ic3);
    
    Index index_child = enzo_block->index().index_child(ic3,min_level_);

    SOLVER_CONTROL(enzo_block,"coarse","fine-1", "23 call prolong_recv");
    enzo::block_array()[index_child].p_solver_mg0_prolong_recv(msg);

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
  TRACE_MG (this,"EnzoBlock::p_solver_mg0_prolong_recv()");

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

  DEBUG_FIELD (enzo_block,ix_,"X prolong_recv");
  // Save message

  // Return if not ready yet
  if (msg != NULL) *pmsg(enzo_block) = msg;
  
  if (! psync_prolong(enzo_block)->next() ) return;
  // Restore saved message then clear
  msg = *pmsg(enzo_block);
  *pmsg(enzo_block) = NULL;

  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "24 prolong_recv");
  
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

  DEBUG_FIELD (enzo_block,ix_,"X update");
  unpack_correction_(enzo_block,msg);
  DEBUG_FIELD (enzo_block,ix_,"X update");

  
  TRACE_LEVEL("EnzoSolverMg0::prolong_recv",enzo_block);
  TRACE_MG (enzo_block,"EnzoSolverMg0::prolong_recv()");
  
  Field field = enzo_block->data()->field();

  enzo_float * X = (enzo_float*) field.values(ix_);
  enzo_float * C = (enzo_float*) field.values(ic_);

  DEBUG_COPY_FIELD("C",ic_);

  DEBUG_FIELD (enzo_block,ix_,"X update");

  for (int i=0; i<mx_*my_*mz_; i++) {
    X[i] += C[i];
  }

  DEBUG_FIELD (enzo_block,ix_,"X update");
  DEBUG_FIELD (enzo_block,ic_,"C update");
  if (index_smooth_post_ >= 0) {

    SOLVER_CONTROL(enzo_block,"coarse","fine-1", "25 call post_smooth");
    call_post_smoother(enzo_block);

  } else {

    SOLVER_CONTROL(enzo_block,"coarse","fine-1", "26 call post_smooth");
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

  solver->post_smooth(this);
  
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
  SOLVER_CONTROL(enzo_block,"coarse+1","fine", "27 post-smooth");

  TRACE_LEVEL("EnzoSolverMg0::post_smooth",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::post_smooth()");

  const int level = enzo_block->level();

  if ( ! is_finest_(enzo_block) ) {

    SOLVER_CONTROL(enzo_block,"coarse","fine-1", "28 call prolong_send");
    prolong_send_ (enzo_block);
  } 

  TRACE_MG(enzo_block,"EnzoSolverMg0::post_smooth() calling end_cycle()");
  SOLVER_CONTROL(enzo_block,"coarse","fine-1", "29 call end_cycle");
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
  SOLVER_CONTROL(enzo_block,"min","max", "30 end_cycle");
  TRACE_LEVEL("EnzoSolverMg0::end_cycle",enzo_block);

  TRACE_MG(enzo_block,"EnzoSolverMg0::end_cycle()");

  ++ (*piter(enzo_block));

  bool is_converged = is_converged_(enzo_block);
  bool is_diverged  = is_diverged_(enzo_block);

  const int iter = *piter(enzo_block);
	    
  const int level = enzo_block->level();

  const bool l_output =
    ( ( enzo_block->index().is_root()) &&
      ( (is_converged) || (is_diverged) ||
	(monitor_iter_ && (iter % monitor_iter_) == 0 )) );

  if (l_output) {
    Solver::monitor_output_(enzo_block,iter,rr0_,0.0,rr_,0.0);
  }

  if (is_converged || is_diverged) {

    // Do an optional final smoothing on the full mesh For use in Dan
    // Reynolds HG algorithm in which Mg0 with no pre- or
    // post-smoothings is used as a preconditioner to BiCgStab
    
    if (index_smooth_last_ >= 0 && (is_finest_(enzo_block)) ) {

      SOLVER_CONTROL(enzo_block,"coarse","fine-1", "31 call last_smooth");
      call_last_smoother(enzo_block);

    } else {

      SOLVER_CONTROL(enzo_block,"coarse","fine-1", "32 call end");
      end (enzo_block);
    }

  } else {

    if ( is_finest_(enzo_block)) {

      SOLVER_CONTROL(enzo_block,"coarse","fine-1", "33 call begin_cycle");
      begin_cycle_ (enzo_block);

    } else {

      int level = enzo_block->level();
      if ( ! (coarse_level_ <= level && level <= max_level_) ) {
	SOLVER_CONTROL(enzo_block,"coarse","fine-1", "34 call coarse_solver");
	DEBUG_FIELD (enzo_block,ix_,"X end_cycle");
	call_coarse_solver(enzo_block);
      }

    }

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

  solver->end(this);
  
#ifdef DEBUG_ENTRY
  if (AFTER_CYCLE(this,CYCLE))
    CkPrintf ("%d %s %p %d DEBUG_ENTRY  exit p_solver_mg0_last_smooth\n",
	      CkMyPe(),name().c_str(),this,index_solver_.back());
#endif
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

FieldMsg * EnzoSolverMg0::pack_residual_(EnzoBlock * enzo_block) throw()
{
  Index index        = enzo_block->index();
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

  Field field = enzo_block->data()->field();

  field_face->face_to_array(field,&narray,&array);

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

  return msg;

}

//----------------------------------------------------------------------

void EnzoSolverMg0::unpack_residual_
(EnzoBlock * enzo_block,FieldMsg * msg) throw()
{
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

  field_face->set_restrict(restrict_);

  Field field = enzo_block->data()->field();
  
  char * a = msg->a;
  field_face->array_to_face(a, field);
  delete field_face;

  delete msg;
}

//----------------------------------------------------------------------

FieldMsg * EnzoSolverMg0::pack_correction_
(EnzoBlock * enzo_block, int ic3[3]) throw()
{
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

  Field field = enzo_block->data()->field();
  field_face->set_prolong(prolong_);

  int narray; 
  char * array;
    
  field_face->face_to_array (field,&narray,&array);

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

#ifdef DEBUG_PROLONG
  if (AFTER_CYCLE(enzo_block,CYCLE))
    CkPrintf ("%s %s %s DEBUG_PROLONG B call prolong_recv()\n",
	      __FILE__,enzo_block->name().c_str(),name_.c_str());
#endif

  return msg;
}

//----------------------------------------------------------------------

void EnzoSolverMg0::unpack_correction_
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
{
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

  field_face->set_prolong(prolong_);

  Field field = enzo_block->data()->field();
  
  field_face->array_to_face (msg->a, field);

  delete field_face;

  // GET DATA
  
  delete msg;

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
  const int iter = *(((EnzoSolverMg0 *)this)->piter(enzo_block));
  return (iter >= iter_max_);
}

//----------------------------------------------------------------------

void EnzoSolverMg0::end(Block * block)
{
  SOLVER_CONTROL(block,"min","max", "35 end");

  TRACE_MG(block,"EnzoSolverMg0::end");
    
  deallocate_temporary_(block);
    
  Solver::end_(block);
}
