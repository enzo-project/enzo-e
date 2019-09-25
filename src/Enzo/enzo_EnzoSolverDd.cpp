
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDd.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-10-01
/// @brief    Implements the EnzoSolverDd class
///
/// @brief [\ref Enzo] Multigrid on the root-level grid using Dd, then
/// BiCgStab in overlapping subdomains defined by root-level Blocks.
/// An optional final Jacobi step can be applied to smooth the solution
/// along subdomain boundaries.
///
/// 0. restrict b to level 0
/// 1. coarse solve: solve Ax=b on level 0
/// 2. Prolong x to child blocks as xc
/// 3. domain solve: solve Ai xi = bi in each root-grid block
///        use xc for boundary conditions and initial guess
/// 4. final smoother: apply final smoother on A x = b (if any)
///
/// 

#include "cello.hpp"
#include "enzo.hpp"

// #define TRACE_DD
// #define TRACE_DD_CYCLE 0

// #define DEBUG_FIELD
// #define DEBUG_COPY_FIELD

#ifdef DEBUG_FIELD

#  define TRACE_FIELD(BLOCK,MSG,ID)					\
  {									\
    Field field = BLOCK->data()->field();				\
    double sum_abs=0.0;							\
    int count = 0;							\
    auto X = (enzo_float*) field.values(ID);				\
    for (int iz=gz_; iz<mz_-gz_; iz++) {				\
      for (int iy=gy_; iy<my_-gy_; iy++) {				\
	for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	  int i = ix + mx_*(iy + my_*iz);				\
	  sum_abs+=std::abs(X[i]);					\
	  count++;							\
	}								\
      }									\
    }									\
    CkPrintf ("%s:%d %s TRACE_FIELD level %d field %d %s %lg count %d\n",		\
	      __FILE__,__LINE__,BLOCK->name().c_str(),BLOCK->level(),ID,MSG,sum_abs,count); \
  }
#else
#  define TRACE_FIELD(BLOCK,MSG,ID) /* ... */
#endif

#ifdef DEBUG_COPY_FIELD
#   define COPY_FIELD(BLOCK,ID,COPY)					\
  {									\
    Field field = BLOCK->data()->field();				\
    enzo_float* X      = (enzo_float*) field.values(ID);		\
    enzo_float* X_bcg  = (enzo_float*) field.values(COPY);		\
    if (X_bcg) for (int i=0; i<mx_*my_*mz_; i++)  X_bcg[i] = X[i];		\
    long double sum_a=0.0,sum_abs=0.0;					\
    for (int iz=gz_; iz<mz_-gz_; iz++) {				\
      for (int iy=gy_; iy<my_-gy_; iy++) {				\
	for (int ix=gx_; ix<mx_-gx_; ix++) {				\
	  int i=ix+mx_*(iy+my_*iz);					\
	  sum_a+=X[i];							\
	  sum_abs+=std::abs(X[i]);					\
	}								\
      }									\
    }									\
    CkPrintf ("%s:%d %s %s COPY_FIELD %d %s shift %Lg %Lg\n"	\
	      ,__FILE__,__LINE__,BLOCK->name().c_str(),name().c_str(),ID,COPY,sum_a, sum_abs); \
  }
#else
#   define COPY_FIELD(BLOCK,ID,COPY) /* ... */
#endif


#ifdef TRACE_DD
#  undef TRACE_DD
#  define TRACE_DD(BLOCK,SOLVER,msg)					\
  if (BLOCK->cycle() >= TRACE_DD_CYCLE) {				\
    CkPrintf ("%d %s %s:%d TRACE_DD %s %s level %d\n",			\
	      CkMyPe(), (BLOCK!=NULL) ? BLOCK->name().c_str():"no block", \
	      __FILE__,__LINE__,					\
	      (SOLVER!=NULL) ? SOLVER->name().c_str():"no solver",msg,	\
	      (BLOCK!=NULL) ? BLOCK->level() : -99);			\
    fflush(stdout);							\
  }
#else
#  define TRACE_DD(BLOCK,SOLVER,msg) /* ... */
#endif


//======================================================================

EnzoSolverDd::EnzoSolverDd
  (std::string name,
   std::string field_x,
   std::string field_b,
   int monitor_iter,
   int restart_cycle,
   int solve_type,
   int min_level,
   int max_level,
   int index_solve_coarse,
   int index_solve_domain,
   int index_solve_smooth,
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
      index_solve_coarse_(index_solve_coarse),
      index_solve_domain_(index_solve_domain),
      index_solve_smooth_(index_solve_smooth),
      restrict_(restrict),
      prolong_(prolong),
      ixc_(-1),
      mx_(0),my_(0),mz_(0),
      gx_(0),gy_(0),gz_(0),
      coarse_level_(coarse_level)
{
  // Initialize temporary fields
  Block * block = NULL;

  ixc_ = cello::field_descr()->insert_temporary();

  /// Initialize default Refresh

  Refresh & refresh = this->refresh_post();
  cello::simulation()->new_refresh_set_name(ir_post_,name);

  refresh.add_field (ix_);
  
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  i_sync_restrict_ = scalar_descr_sync->new_value(name + ":restrict");
  i_sync_prolong_  = scalar_descr_sync->new_value(name + ":prolong");

  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  i_msg_ = scalar_descr_void->new_value(name + ":msg");
}

//----------------------------------------------------------------------

void EnzoSolverDd::apply ( std::shared_ptr<Matrix> A, Block * block) throw()
{
  TRACE_DD(block,this,"apply");

  Solver::begin_(block);

  Field field = block->data()->field();

  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  TRACE_FIELD(block,"apply B",ib_);

  A_ = A;

  allocate_temporary_(block);

  // Check that component solvers are of the correct type
  ASSERT2("EnzoSolverDd::apply()",
	  "Coarse solver %s type %s != solve_level",
	  cello::solver(index_solve_coarse_)->name().c_str(),
	  solve_string[cello::solver(index_solve_coarse_)->solve_type()],
	  cello::solver(index_solve_coarse_)->solve_type() == solve_level);
  
  Sync * sync_restrict = psync_restrict(block);

  sync_restrict->reset();
  sync_restrict->set_stop(cello::num_children());
  
  Sync * sync_prolong = psync_prolong(block);

  sync_prolong->reset();
  sync_prolong->set_stop(2); // self and parent

  int level = block->level();

  const int m = mx_*my_*mz_;
  
  if ( ! block->is_leaf() ) {
    std::fill_n ((enzo_float*) field.values(ib_), m, 0.0);
  }

  std::fill_n ((enzo_float*) field.values(ix_),  m, 0.0);
  std::fill_n ((enzo_float*) field.values(ixc_), m, 0.0);
	  
  if (block->is_leaf()) {

    COPY_FIELD(block,ib_,"B1_dd");
    begin_solve(enzo::block(block));

  } else {

    const bool in_range = (coarse_level_ <= level && level <= max_level_);
    
    if ( ! in_range ) {
    
      call_coarse_solver(enzo::block(block));
    
    }
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::begin_solve(EnzoBlock * enzo_block) throw()
{
  TRACE_DD(enzo_block,this,"begin_solve");
  TRACE_FIELD(enzo_block,"begin_solve B",ib_);

  if (enzo_block->level() == coarse_level_) {

    call_coarse_solver(enzo_block);

  } else {

    restrict (enzo_block);

  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::restrict(EnzoBlock * enzo_block) throw()
{
  TRACE_DD(enzo_block,this,"restrict");

  restrict_send(enzo_block);

  call_coarse_solver(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverDd::restrict_send(EnzoBlock * enzo_block) throw()
{
  TRACE_DD(enzo_block,this,"restrict_send");

  // Pack field
  Index index = enzo_block->index();
  int level   = index.level();  
  int ic3[3];
  index.child(level,&ic3[0],&ic3[1],&ic3[2],min_level_);

  TRACE_FIELD(enzo_block,"restrict_send B",ib_);
  
  FieldMsg * msg = pack_field_(enzo_block,ib_,refresh_coarse,ic3);

  // Send packed field to parent
  Index index_parent = enzo_block->index().index_parent(min_level_);
  enzo::block_array()[index_parent].p_solver_dd_restrict_recv(msg);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_restrict_recv(FieldMsg * msg)
{
  TRACE_DD(this,solver(),"p_solver_dd_restrict_recv");

  static_cast<EnzoSolverDd*> (solver())->restrict_recv(this,msg);
}

//----------------------------------------------------------------------

void EnzoSolverDd::restrict_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
{
  // Unpack "B" vector data from children

  TRACE_DD(enzo_block,this,"restrict_recv");
  unpack_field_(enzo_block,msg,ib_,refresh_coarse);
  TRACE_FIELD(enzo_block,"restrict_recv B",ib_);
  
  if (psync_restrict(enzo_block)->next())  {
    begin_solve(enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::call_coarse_solver(EnzoBlock * enzo_block) throw()
{
  TRACE_DD(enzo_block,this,"call_coarse_solver");
  Solver * solve_coarse = cello::solver(index_solve_coarse_);

  solve_coarse->set_sync_id (enzo_sync_id_solver_dd_coarse);
  solve_coarse->set_callback(CkIndex_EnzoBlock::p_solver_dd_solve_coarse());
  
  COPY_FIELD(enzo_block,ib_,"B2_dd");
  TRACE_FIELD(enzo_block,"coarse_solve B",ib_);
  solve_coarse->set_field_x (ixc_);
  solve_coarse->set_field_b (ib_);

  solve_coarse->apply(A_,enzo_block);
}


//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_solve_coarse()
{
  TRACE_DD(this,solver(),"p_solver_dd_solve_coarse");

  CkCallback callback(CkIndex_EnzoBlock::r_solver_dd_barrier(NULL), 
		      enzo::block_array());
  contribute(callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_dd_barrier(CkReductionMsg * msg)
{
  static_cast<EnzoSolverDd*> (solver())->prolong(this);
  delete msg;
}

//----------------------------------------------------------------------

void EnzoSolverDd::prolong(EnzoBlock * enzo_block) throw()
{

  if (is_finest_(enzo_block)) {
    copy_xc_to_x_(enzo_block);
  }

  /// Prolong solution to next-finer level

  const int level = enzo_block->level();

  if (level == coarse_level_) {

    prolong_send_ (enzo_block);
      
  }
  
  if (coarse_level_ < level && level <= max_level_) {
    TRACE_DD(enzo_block,this,"call_solver_dd_prolong_recv_A");
    enzo_block->solver_dd_prolong_recv(NULL);
  } else {
    TRACE_DD(enzo_block,this,"call_call_domain_solver_A");
    call_domain_solver (enzo_block);
  }
}

//----------------------------------------------------------------------

void EnzoSolverDd::prolong_send_(EnzoBlock * enzo_block) throw()
{
  if ( ! is_finest_(enzo_block) ) {
    TRACE_DD(enzo_block,this,"prolong_send");

    ItChild it_child(cello::rank());
    int ic3[3];

    while (it_child.next(ic3)) {

      TRACE_FIELD(enzo_block,"prolong_send XC",ixc_);
      FieldMsg * msg = pack_field_(enzo_block,ixc_,refresh_fine,ic3);
    
      Index index_child = enzo_block->index().index_child(ic3,min_level_);

      TRACE_DD(enzo_block,this,"call_solver_dd_prolong_recv_B");
      enzo::block_array()[index_child].p_solver_dd_prolong_recv(msg);

    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_prolong_recv(FieldMsg * msg)
{  solver_dd_prolong_recv(msg); }

void EnzoBlock::solver_dd_prolong_recv(FieldMsg * msg)
{
  TRACE_DD(this,solver(),"p_solver_dd_prolong_recv");
  static_cast<EnzoSolverDd*> (solver())->prolong_recv(this,msg);
}

//----------------------------------------------------------------------

void EnzoSolverDd::prolong_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
{
  
  TRACE_DD(enzo_block,this,"prolong_recv_1");
  // Save message
  if (msg != NULL) *pmsg(enzo_block) = msg;
  TRACE_DD(enzo_block,this,"prolong_recv_2");

  // Return if not ready yet
  if (! psync_prolong(enzo_block)->next() ) return;
  TRACE_DD(enzo_block,this,"prolong_recv_3");

  // Restore saved message then clear
  msg = *pmsg(enzo_block);
  *pmsg(enzo_block) = NULL;

  unpack_field_(enzo_block,msg,ixc_,refresh_fine);
  TRACE_FIELD(enzo_block,"prolong_recv XC",ixc_);

  // copy X = XC
  // copy X_copy = XC (using Solver::reuse_solution_(cycle) )
  copy_xc_to_x_(enzo_block);
  
  TRACE_FIELD(enzo_block,"prolong_recv X",ix_);

  COPY_FIELD(enzo_block,ixc_,"XC2_dd");
  TRACE_DD(enzo_block,this,"call_call_domain_solver_B");
  call_domain_solver(enzo_block);

  prolong_send_ (enzo_block);

  TRACE_DD(enzo_block,this,"call_call_domain_solver_C");
}

//----------------------------------------------------------------------

void EnzoSolverDd::copy_xc_to_x_(EnzoBlock * enzo_block) throw()
{

  const int m = mx_*my_*mz_;

  Field field = enzo_block->data()->field();

  std::copy_n((enzo_float *) field.values(ixc_),m,
	      (enzo_float *) field.values(ix_));
  std::copy_n((enzo_float *) field.values(ixc_),m,
	      (enzo_float *) field.values("X_copy"));

}

//----------------------------------------------------------------------
  
void EnzoSolverDd::call_domain_solver(EnzoBlock * enzo_block) throw()
{
  TRACE_DD(enzo_block,this,"call_domain_solver");
  Solver * solve_domain = cello::solver(index_solve_domain_);

  solve_domain->set_sync_id (enzo_sync_id_solver_dd_domain);
  solve_domain->set_callback(CkIndex_EnzoBlock::p_solver_dd_solve_domain());

  solve_domain->set_field_x (ix_);
  solve_domain->set_field_b (ib_);

  COPY_FIELD(enzo_block,ib_,"B3_dd");
  COPY_FIELD(enzo_block,ix_,"X1_dd");
  solve_domain->apply(A_,enzo_block);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_solve_domain()
{
  TRACE_DD(this,solver(),"p_solver_dd_solve_domain");
  static_cast<EnzoSolverDd*> (solver())->continue_after_domain_solve(this);
}

//----------------------------------------------------------------------

void EnzoSolverDd::continue_after_domain_solve(EnzoBlock * enzo_block) throw()
{
  TRACE_FIELD(enzo_block,"domain_solve XC",ixc_);
  TRACE_FIELD(enzo_block,"domain_solve X",ix_);
  TRACE_DD(enzo_block,this,"continue_after_domain_solve");
  COPY_FIELD(enzo_block,ix_,"X2_dd");
  CkCallback callback(CkIndex_EnzoBlock::r_solver_dd_end(NULL), 
		      enzo::block_array());
  enzo_block->contribute(callback);

}

//----------------------------------------------------------------------

void EnzoBlock::r_solver_dd_end(CkReductionMsg * msg)
{
  static_cast<EnzoSolverDd*> (solver())->call_last_smoother(this);
  delete msg;
}

//----------------------------------------------------------------------

void EnzoSolverDd::call_last_smoother(EnzoBlock * enzo_block) throw()
{
  TRACE_DD(enzo_block,this,"call_last_smoother");
  Solver * smooth_last = cello::solver(index_solve_smooth_);

  smooth_last->set_sync_id (enzo_sync_id_solver_dd_smooth);
  smooth_last->set_callback(CkIndex_EnzoBlock::p_solver_dd_last_smooth());

  smooth_last->set_field_x(ix_);
  smooth_last->set_field_b(ib_);
  COPY_FIELD(enzo_block,ib_,"B4_dd");

  TRACE_FIELD(enzo_block,"end",ix_);
  smooth_last->apply(A_,enzo_block);

}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_dd_last_smooth()
{
  TRACE_DD(this,solver(),"p_solver_dd_last_smooth");
  static_cast<EnzoSolverDd*> (solver())->continue_after_last_smooth(this);
}

//----------------------------------------------------------------------

void EnzoSolverDd::continue_after_last_smooth(EnzoBlock * enzo_block) throw()
{
  TRACE_FIELD(enzo_block,"last_smooth XC",ixc_);
  TRACE_FIELD(enzo_block,"last_smooth X",ix_);
  TRACE_DD(enzo_block,this,"continue_after_last_smooth");
  COPY_FIELD(enzo_block,ix_,"X3_dd");

  end(enzo_block);
}

//----------------------------------------------------------------------

void EnzoSolverDd::end (Block* block) throw ()
{
  TRACE_FIELD(block,"end",ix_);
  TRACE_DD(block,this,"end");
  deallocate_temporary_(block);
  
  Solver::end_(block);
}

//======================================================================

FieldMsg * EnzoSolverDd::pack_field_(EnzoBlock * enzo_block,
				     int index_field,
				     int refresh_type,
				     int * ic3)
{
  int  if3[3] = {0,0,0};
  bool lg3[3];
  for (int i=0; i<3; i++) lg3[i] = (refresh_type == refresh_fine);

  Refresh * refresh = new Refresh;
  refresh->add_field(index_field);

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, lg3, refresh_type, refresh, true);

  if (refresh_type == refresh_coarse)
    field_face->set_restrict(restrict_);
  if (refresh_type == refresh_fine)
    field_face->set_prolong(prolong_);
  
  Field field = enzo_block->data()->field();
  int narray; 
  char * array;
  field_face->face_to_array(field,&narray,&array);

  delete field_face;

  FieldMsg * msg  = new (narray) FieldMsg;
 
  msg->n = narray;
  memcpy (msg->a, array, narray);
  delete [] array;
  
  msg->ic3[0] = ic3[0];
  msg->ic3[1] = ic3[1];
  msg->ic3[2] = ic3[2];

  return msg;

}

//----------------------------------------------------------------------

void EnzoSolverDd::unpack_field_
(EnzoBlock * enzo_block,
 FieldMsg * msg,
 int index_field,
 int refresh_type)
{
  int if3[3] = {0,0,0};
  bool lg3[3];
  for (int i=0; i<3; i++) lg3[i] = (refresh_type == refresh_fine);
  Refresh * refresh = new Refresh;
  refresh->add_field(index_field);

  int * ic3 = msg->ic3;

  FieldFace * field_face = enzo_block->create_face 
    (if3, ic3, lg3, refresh_type, refresh, true);

  if (refresh_type == refresh_coarse)
    field_face->set_restrict(restrict_);
  if (refresh_type == refresh_fine)
    field_face->set_prolong(prolong_);

  Field field = enzo_block->data()->field();
  
  char * a = msg->a;
  field_face->array_to_face(a, field);
  delete field_face;

  delete msg;
}

