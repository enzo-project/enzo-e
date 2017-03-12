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
///      if (sync.next())
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

// #define DEBUG_SOLVER_MG0

// #define DEBUG_TRACE_LEVEL

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
    double xx=0;							\
    for (int iz=0; iz<mz_; iz++) {					\
      for (int iy=0; iy<my_; iy++) {					\
	for (int ix=0; ix<mx_; ix++) {					\
	  int i = ix + mx_*(iy + my_*iz);				\
	  xx+=X[i]*X[i];						\
	}								\
      }									\
    }									\
    CkPrintf ("%s:%d %s DEBUG_SOLVER ||%s|| = %g\n",			\
	      __FILE__,__LINE__,enzo_block->name().c_str(),NAME,xx);	\
  }

#else
#   define TRACE_MG(block,msg) /*  ... */

#   define DEBUG_FIELD(IX,NAME) /* ... */
#endif

#ifdef DEBUG_SOLVER_MG0
#define LOCAL_SUM(M,it)							\
  {									\
    T * X = (T*) enzo_block->data()->field().values(it);		\
    double min=std::numeric_limits<double>::max();			\
    double max=-std::numeric_limits<double>::max();			\
    double sum =0.0;							\
    const int i0 = gx_ + mx_*(gy_ + my_*gz_);				\
									\
    for (int iz=0; iz<nz_; iz++) {					\
      for (int iy=0; iy<ny_; iy++) {					\
	for (int ix=0; ix<nx_; ix++) {					\
	  int i = i0 + ix + mx_*(iy + my_*iz);				\
	  min = std::min(min,double(X[i]));				\
	  max = std::max(max,double(X[i]));				\
	  sum += X[i];							\
	}								\
      }									\
    }									\
    CkPrintf ("LOCAL_SUM %d iter %d level %d %s block %s"		\
	      "min %20.15e sum %20.15e max %20.15e\n",			\
	      __LINE__,enzo_block->mg_iter(), enzo_block->level(),M,	\
	      enzo_block->name().c_str(), min,sum,max);			\
  }
#else
#define LOCAL_SUM(M,it) /* EMPTY */
#endif  

//======================================================================

extern CkReduction::reducerType sum_long_double_2_type;

//======================================================================

EnzoSolverMg0::EnzoSolverMg0 
(const FieldDescr * field_descr, 
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
    A_(new EnzoMatrixLaplace),
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
  if (min_level_ > 0) {
    ERROR1 ("EnzoSolverMg0::EnzoSolverMg0()",
	    "Solver Mg0 requires min_level = %d to be <= 0",
	    min_level_);
  }

  /// Initialize default Refresh

  add_refresh(4,0,neighbor_level,sync_barrier);
  refresh(0)->add_all_fields(field_descr->field_count());

  ir_ = field_descr->field_id("R");
  ic_ = field_descr->field_id("C");

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

  TRACE_LEVEL("EnzoSolverMg0::apply",block);

  Field field = block->data()->field();

  field.size           (&nx_,&ny_,&nz_);
  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  // Initialize child counter for restrict synchronization
  enzo_block->mg_sync_set_stop(NUM_CHILDREN(enzo_block->rank()));
  enzo_block->mg_sync_reset();

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

    CkCallback callback(CkIndex_EnzoBlock::p_solver_mg0_shift_b<T>(NULL), 
			enzo_block->proxy_array());

    enzo_block->contribute(2*sizeof(long double), &reduce, 
			   sum_long_double_2_type, callback);

  } else {
    begin_solve<T>(enzo_block);
  }

}

//----------------------------------------------------------------------

template<class T>
void EnzoBlock::p_solver_mg0_shift_b(CkReductionMsg* msg)
{

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
  solver->begin_solve<T>(this);

}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::begin_solve(EnzoBlock * enzo_block) throw()
{
  TRACE_LEVEL("EnzoSolverMg0::begin_solve",enzo_block);
  // start the MG V-cycle with the max level blocks

  if (enzo_block->is_leaf() || (enzo_block->level() == max_level_)) {

    if (A_->is_singular()) {

      Field field = enzo_block->data()->field();
      T* B  = (T*) field.values(ib_);
      T shift = -bs_ / bc_;
#ifdef DEBUG_SOLVER_MG0      
      CkPrintf ("%d %s DEBUG_SOLVER_MG0 bs %lf bc %lf\n",
		__LINE__,enzo_block->name().c_str(),double(bs_),double(bc_));
#endif      

      LOCAL_SUM("B shift 0",ib_);
      for (int iz=0; iz<mz_; iz++) {
	for (int iy=0; iy<my_; iy++) {
	  for (int ix=0; ix<mx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    B[i] += shift;
	  }
	}
      }

      LOCAL_SUM("B shift 1",ib_);
    }

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

  if (is_converged_(enzo_block)) {

    TRACE_MG(enzo_block,"converged");

    // CkCallback(callback_,
    // 	       CkArrayIndexIndex(enzo_block->index()),
    // 	       enzo_block->proxy_array()).send();

  } else if (level == min_level_) {
    // TODO REFRESH X
    TRACE_MG(enzo_block,"calling coarse solve");

    Field field = enzo_block->data()->field();
    T * X = (T*) field.values(ix_);
    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  X[i] = 0.0;
	}
      }
    }

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
      for (int iz=0; iz<mz_; iz++) {
	for (int iy=0; iy<my_; iy++) {
	  for (int ix=0; ix<mx_; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    X[i] = 0;
	  }
	}
      }
    }

#ifdef DEBUG_SOLVER_REFRESH    
    CkPrintf ("DEBUG_SOLVER_MG refresh sync_face %d\n",refresh.sync_id());
#endif

    Data * data = enzo_block->data();
    Field field = data->field();

    Simulation * simulation = proxy_simulation.ckLocalBranch();

    if (index_smooth_pre_ >= 0) {

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
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_mg0_pre_smooth()
/// [*]
{
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
  
  solve_coarse->apply(A_,ix_,ib_,enzo_block);

  
}

//----------------------------------------------------------------------

// template <class T>
// void EnzoBlock::p_solver_mg0_restrict_send(CkReductionMsg * msg)
// /// [*]
// {
//   TRACE_MG(this,"EnzoBlock::p_solver_mg0_restrict_send()");

//   EnzoSolverMg0 * solver = 
//     static_cast<EnzoSolverMg0*> (this->solver());

//   // GET DATA
  
//   delete msg;

//   solver -> restrict_send<T>(this);
// }

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
  int level          = index.level();
  
  // copy face data to FieldFace

  // Pack and send "R" to parent

  int ic3[3];
  index.child(level,&ic3[0],&ic3[1],&ic3[2],min_level_);

  // <COMMON CODE> in restrict_send_() and prolong_send_()
  
  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  std::vector<int> field_list;
  field_list.push_back(ir_);

  // copy data from EnzoBlock to array via FieldFace

  FieldFace * field_face = enzo_block->create_face
    (if3, ic3, lg3, refresh_coarse, field_list);

  field_face->set_restrict(restrict_);
  
  int narray; 
  char * array;

  field_face->face_to_array(enzo_block->data()->field(),&narray,&array);

  delete field_face;

  // Create a FieldMsg for sending data to parent
  // (note: charm messages not deleted on send; are deleted on receive)

  FieldMsg * field_message  = new (narray) FieldMsg;
 
  /// WARNING: double copy

  // Copy FieldFace data to field_message

  field_message->n = narray;
  memcpy (field_message->a, array, narray);
  field_message->ic3[0] = ic3[0];
  field_message->ic3[1] = ic3[1];
  field_message->ic3[2] = ic3[2];

  //  </COMMON CODE>

  DEBUG_FIELD(ir_,"R");
  
  CkCallback (CkIndex_EnzoBlock::p_solver_mg0_restrict_recv<T>(NULL), 
	      CkArrayIndexIndex(index_parent),
	      enzo_block->proxy_array()).send(field_message);

  delete [] array;

}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_mg0_restrict_recv(FieldMsg * field_message)
/// [*]
{

  TRACE_MG(this,"EnzoBlock::restrict_recv()");

  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  solver -> restrict_recv<T>(this,field_message);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::restrict_recv(EnzoBlock * enzo_block, FieldMsg * field_message) throw()
/// 
/// [*] restrict recv
///
///      [ unpack B ]
///      if (sync.next())
///          begin_cycle()
{

  // Unpack "B" vector data from children

  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  std::vector<int> field_list;
  field_list.push_back(ib_);

  // copy data from field_message to this EnzoBlock

  int * ic3 = field_message->ic3;

  FieldFace * field_face = enzo_block->create_face 
    (if3, ic3, lg3, refresh_coarse, field_list);

  field_face->set_restrict(restrict());

  char * a = field_message->a;
  field_face->array_to_face(a, enzo_block->data()->field());
  delete field_face;

  delete field_message;

  // continue with EnzoSolverMg0


  TRACE_LEVEL("EnzoSolverMg0::restrict_recv",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_recv()");

  DEBUG_FIELD(ib_,"B");
  
  if (enzo_block->mg_sync_next()) {
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
  DEBUG_FIELD(ix_,"X");
  while (it_child.next(ic3)) {

    Index index_child = enzo_block->index().index_child(ic3,min_level_);

    // Pack and send "X" to children

    // <COMMON CODE> in restrict_send_() and prolong_send_()

    int if3[3] = {0,0,0};
    bool lg3[3] = {false,false,false};
    std::vector<int> field_list;
    field_list.push_back(ix_);

    // copy data from EnzoBlock to array via FieldFace

    FieldFace * field_face = enzo_block->create_face
      (if3, ic3, lg3, refresh_fine, field_list);

    field_face->set_prolong(prolong_);

    int narray; 
    char * array;
    
    field_face->face_to_array (enzo_block->data()->field(),&narray,&array);

    delete field_face;

    // Create a FieldMsg for sending data to parent
    // (note: charm messages not deleted on send; are deleted on receive)
    
    FieldMsg * field_message  = new (narray) FieldMsg;

    /// WARNING: double copy

    // Copy FieldFace data to field_message

    field_message->n = narray;
    memcpy (field_message->a, array, narray);
    field_message->ic3[0] = ic3[0];
    field_message->ic3[1] = ic3[1];
    field_message->ic3[2] = ic3[2];

    //  </COMMON CODE>

    CkCallback (CkIndex_EnzoBlock::p_solver_mg0_prolong_recv<T>(NULL), 
		CkArrayIndexIndex(index_child),
		enzo_block->proxy_array()).send(field_message);

    delete [] array;

  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg * field_message)
/// [*]
{

  TRACE_MG (this,"EnzoBlock::p_solver_mg0_prolong_recv()");
  
  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  solver -> prolong_recv<T>(this,field_message);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::prolong_recv(EnzoBlock * enzo_block, FieldMsg * field_message) throw()
/// 
/// [ ] prolong recv
///
///      [ unpack C ]
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{

  // Unpack "C" vector data from children

  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  std::vector<int> field_list;
  field_list.push_back(ic_);

  // copy data from field_message to this EnzoBlock

  FieldFace * field_face = enzo_block->create_face 
    (if3, field_message->ic3, lg3, refresh_fine, field_list);

  field_face->set_prolong(prolong());

  field_face->array_to_face (field_message->a, enzo_block->data()->field());

  delete field_face;

  // GET DATA
  
  delete field_message;

  DEBUG_FIELD(ic_,"C");

  TRACE_LEVEL("EnzoSolverMg0::prolong_recv",enzo_block);
  TRACE_MG (enzo_block,"EnzoSolverMg0::prolong_recv()");
  
  Data * data = enzo_block->data();
  Field field = data->field();

  T * X = (T*) field.values(ix_);
  T * C = (T*) field.values(ic_);

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

  if (is_converged_(enzo_block)) {

    Solver::end_(enzo_block);
    CkCallback(callback_,
    	       CkArrayIndexIndex(enzo_block->index()),
    	       enzo_block->proxy_array()).send();

  }

  if (enzo_block->is_leaf() || (enzo_block->level() == max_level_)) {

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

//======================================================================
