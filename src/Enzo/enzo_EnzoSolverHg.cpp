// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverHg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoSolverHg class
///
/// Multigrid-like preconditioner for adaptive meshes developed by
/// Prof. Daniel Reynolds.  No smoothings are performed, just a
/// restriction of the fine-grid rhs to the root grid, a root-level
/// solve, then projection of the root-level solution back to the
/// fine-level grid.  Could include smoothings but without refresh.
///
///======================================================================
///
///  "Coarse" view of HG preconditioner
///
///   @code
///
///   $HG(A_h,X_h,B_h)$
///
///       if (level == min_level) then
///          solve_coarse()     solve $A_h X_h = B_h$
///       else
///          p_restrict ()      $B_H = I_h^H B_h$
///          MG()               $HG(A_H, X_H, B_H)$
///          p_prolong ()       $X_h = I_H^h X_H$
///
///  @endcode
///
///----------------------------------------------------------------------
///
///  "Fine" view of HG multigrid solver
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

// #define DEBUG_SOLVER_HG

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
  

#ifdef DEBUG_SOLVER_HG
#   define TRACE_MG(block,msg)						\
  CkPrintf ("%d %s TRACE_MG %s\n",					\
	    CkMyPe(),(block != NULL) ? block->name().c_str() : "root",msg); \
  fflush(stdout);

#define DEBUG_PRINT_XX							\
  {									\
    Data * data = enzo_block->data();					\
    Field field = data->field();					\
    T * X = (T*) field.values(ix_);					\
    double xx=0;							\
    for (int iz=0; iz<mz_; iz++) {					\
      for (int iy=0; iy<my_; iy++) {					\
	for (int ix=0; ix<mx_; ix++) {					\
	  int i = ix + mx_*(iy + my_*iz);				\
	  xx+=X[i]*X[i];						\
	}								\
      }									\
    }									\
    CkPrintf ("%d DEBUG_SOLVER %s ix=%d xx=%g\n",__LINE__,enzo_block->name().c_str(),ix_); \
  }

#else
#   define TRACE_MG(block,msg) /*  ... */

#   define DEBUG_PRINT_XX /* ... */
#endif

#ifdef DEBUG_SOLVER_HG
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

EnzoSolverHg::EnzoSolverHg 
(const FieldDescr * field_descr, 
 int monitor_iter,
 int rank,
 int iter_max,
 int index_solve_coarse,
 Restrict * restrict,
 Prolong * prolong,
 int min_level,
 int max_level) 
  : Solver(), 
    A_(new EnzoMatrixLaplace),
    index_solve_coarse_(index_solve_coarse),
    restrict_(restrict),
    prolong_(prolong),
    rank_(rank),
    iter_max_(iter_max),
    monitor_iter_(monitor_iter),
    ib_(0), ic_(0), ir_(0), ix_(0),
    min_level_(min_level),
    max_level_(max_level),
    mx_(0),my_(0),mz_(0),
    nx_(0),ny_(0),nz_(0),
    gx_(0),gy_(0),gz_(0),
    bs_(0), bc_(0)
    // [*]
{
  if (min_level_ > 0) {
    ERROR1 ("EnzoSolverHg::EnzoSolverHg()",
	    "Solver Hg requires min_level = %d to be <= 0",
	    min_level_);
  }

  /// Initialize default Refresh

  add_refresh(4,0,neighbor_level,sync_barrier);
  refresh(0)->add_all_fields(field_descr->field_count());

  ir_ = field_descr->field_id("R");
  ic_ = field_descr->field_id("C");
  id_ = field_descr->field_id("D");

}

//----------------------------------------------------------------------

EnzoSolverHg::~EnzoSolverHg () throw()
// [*]
{
  delete prolong_;
  delete restrict_;

  prolong_ = NULL;
  restrict_ = NULL;
}


//----------------------------------------------------------------------

void EnzoSolverHg::apply
( Matrix * A, int ix, int ib, Block * block) throw()
{
  TRACE_MG(block,"EnzoSolverHg::apply()");

  block->push_solver(index_);

  A_ = A;
  ix_ = ix;
  ib_ = ib;

  TRACE_LEVEL("EnzoSolverHg::apply",block);

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
    ERROR1("EnzoSolverHg()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::enter_solver_ (EnzoBlock * enzo_block) throw()
/// [*]
///
///     iter = 0
///     initialize X,B,R,C
///     if (level == max_level) 
///        begin_cycle()
{
  TRACE_LEVEL("EnzoSolverHg::enter_solver",enzo_block);

  TRACE_MG(enzo_block,"EnzoSolverHg::enter_solver()");

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

#ifdef DEBUG_SOLVER_HG    
    CkPrintf ("%d %s DEBUG_SOLVER_HG bs %lf bc %lf\n",
	      __LINE__,enzo_block->name().c_str(),double(reduce[0]),double(reduce[1]));
#endif    

    /// initiate callback for r_solver_bicgstab_start_1 and
    /// contribute to sum and count

    // CkPrintf ("%s %d p_solver_hg_shift_b\n",
    // 	      enzo_block->name().c_str(), enzo_block->level());
    // fflush(stdout);
    CkCallback callback(CkIndex_EnzoBlock::p_solver_hg_shift_b<T>(NULL), 
			enzo_block->proxy_array());
    
    enzo_block->contribute(2*sizeof(long double), &reduce, 
			   sum_long_double_2_type, callback);

  } else {
    begin_solve<T>(enzo_block);
  }

}

//----------------------------------------------------------------------

template<class T>
void EnzoBlock::p_solver_hg_shift_b(CkReductionMsg* msg)
{

  /// EnzoBlock accumulates global contributions to SUM(B) and COUNT(B)
  EnzoSolverHg* solver = 
    static_cast<EnzoSolverHg*> (this->solver());
  
  long double* data = (long double*) msg->getData();

  solver->set_bs( data[0] );
  solver->set_bc( data[1] );

#ifdef DEBUG_SOLVER_HG  
  CkPrintf ("%d %s DEBUG_SOLVER_HG bs %lf bc %lf\n",
	    __LINE__,name().c_str(),double(data[0]),double(data[1]));
#endif  

  delete msg;

  /// Start solver
  solver->begin_solve<T>(this);

}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::begin_solve(EnzoBlock * enzo_block) throw()
{
  TRACE_LEVEL("EnzoSolverHg::begin_solve",enzo_block);
  // start the MG V-cycle with the max level blocks

  if (enzo_block->level() == max_level_) {

    if (A_->is_singular()) {

      Field field = enzo_block->data()->field();
      T* B  = (T*) field.values(ib_);
      T shift = -bs_ / bc_;
#ifdef DEBUG_SOLVER_HG      
      CkPrintf ("%d %s DEBUG_SOLVER_HG bs %lf bc %lf\n",
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
void EnzoSolverHg::begin_cycle_(EnzoBlock * enzo_block) throw()
/// [*]
///
///     if (converged()) exit()
///     if (level == min_level) then
///        coarse_solve(A,X,B)
///     else
///        callback = p_pre_smooth()
///        call refresh (X,"level")
{
  TRACE_LEVEL("EnzoSolverHg::begin_cycle",enzo_block);
  const int level = enzo_block->level();

  //  LOCAL_SUM("0 X",ix_);
  //  LOCAL_SUM("0 R",ir_);
  //  LOCAL_SUM("0 B",ib_);
  //  LOCAL_SUM("0 C",ic_);
  
  TRACE_MG(enzo_block,"EnzoSolverHg::begin_cycle()");

  if (is_converged_(enzo_block)) {

    TRACE_MG(enzo_block,"converged");
    // implicitly exit solver

  } else if (level == min_level_) {
    // TODO REFRESH X
    TRACE_MG(enzo_block,"calling coarse solve");

    const int field_count = enzo_block->data()->field().field_count();

    Refresh refresh (4,0,neighbor_level, sync_face, 4);
    refresh.add_all_fields(field_count);

    enzo_block->refresh_enter
      (CkIndex_EnzoBlock::p_solver_hg_solve_coarse(),&refresh);

    //    solve_coarse_<T>(enzo_block);
    //    LOCAL_SUM("4 X",ix_);

  } else {

    TRACE_MG(enzo_block,"calling smoother");

    const int field_count = enzo_block->data()->field().field_count();

    Refresh refresh (4,0,neighbor_level, sync_face, 6);
    refresh.add_all_fields(field_count);

    enzo_block->refresh_enter
      (CkIndex_EnzoBlock::p_solver_hg_pre_smooth(),&refresh);

  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_hg_solve_coarse()
/// [*]
{
  TRACE_MG(this,"EnzoBlock::solver_hg_coarse_solve()");
  
  EnzoSolverHg * solver = 
    static_cast<EnzoSolverHg*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  if      (precision == precision_single)    
    solver->solve_coarse<float>(enzo_block);
  else if (precision == precision_double)    
    solver->solve_coarse<double>(enzo_block);
  else if (precision == precision_quadruple) 
    solver->solve_coarse<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverBiCgStab()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

void EnzoBlock::p_solver_hg_pre_smooth()
/// [*]
{
  TRACE_MG(this,"EnzoBlock::solver_hg_pre_smooth()");
  
  EnzoSolverHg * solver = 
    static_cast<EnzoSolverHg*> (this->solver());

  EnzoBlock* enzo_block = static_cast<EnzoBlock*> (this);
  Field field = data()->field();
  int precision = field.precision(field.field_id("density")); // assuming 
  if      (precision == precision_single)    
    solver->pre_smooth<float>(enzo_block);
  else if (precision == precision_double)    
    solver->pre_smooth<double>(enzo_block);
  else if (precision == precision_quadruple) 
    solver->pre_smooth<long double>(enzo_block);
  else 
    ERROR1("EnzoSolverBiCgStab()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::pre_smooth(EnzoBlock * enzo_block) throw()
/// [*]
///
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  TRACE_LEVEL("EnzoSolverHg::pre_smooth",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverHg::pre_smooth()");

  Data * data = enzo_block->data();
  Field field = data->field();

  //  smooth_pre_->apply(A_,ix_,ib_,enzo_block);

  restrict_send<T>(enzo_block);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_hg_restrict_send(CkReductionMsg * msg)
/// [*]
{
  TRACE_MG(this,"EnzoBlock::p_solver_hg_restrict_send()");

  EnzoSolverHg * solver = 
    static_cast<EnzoSolverHg*> (this->solver());

  // GET DATA
  
  delete msg;

  solver -> restrict_send<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::restrict_send(EnzoBlock * enzo_block) throw()
/// 
/// [*] restrict send
///
///      A.residual(R,B,X)
///      pack R
///      index_parent.p_restrict_recv(R)
{
  TRACE_LEVEL("EnzoSolverHg::restrict_send",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverHg::restrict_send()");

  Data * data = enzo_block->data();
  Field field = data->field();

  A_->residual(ir_, ib_, ix_, enzo_block);

  //  LOCAL_SUM("2 R",ir_);

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

  // CkPrintf ("%s p_solver_hg_restrict_recv\n",
  // 	    enzo_block->name().c_str(), enzo_block->level());
  // fflush(stdout);
  CkCallback (CkIndex_EnzoBlock::p_solver_hg_restrict_recv<T>(NULL), 
	      CkArrayIndexIndex(index_parent),
	      enzo_block->proxy_array()).send(field_message);

  delete [] array;

}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_hg_restrict_recv(FieldMsg * field_message)
/// [*]
{

  TRACE_MG(this,"EnzoBlock::restrict_recv()");

  // Unpack "B" vector data from children

  // Face (0,0,0) is the entire block
  int if3[3] = {0,0,0};
  // exclude ghost zones
  bool lg3[3] = {false,false,false};
  // receive in vector "B"
  std::vector<int> field_list;
  field_list.push_back(data()->field().field_id("B"));

  // copy data from field_message to this EnzoBlock

  EnzoSolverHg * solver = 
    static_cast<EnzoSolverHg*> (this->solver());

  int * ic3 = field_message->ic3;

  FieldFace * field_face = create_face 
    (if3, ic3, lg3, refresh_coarse, field_list);

  field_face->set_restrict(solver->restrict());

  char * a = field_message->a;
  field_face->array_to_face(a, data()->field());
  delete field_face;

  delete field_message;

  // continue with EnzoSolverHg

  solver -> restrict_recv<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::restrict_recv(EnzoBlock * enzo_block) throw()
/// 
/// [*] restrict recv
///
///      [ unpack B ]
///      if (sync.next())
///          begin_cycle()
{
  TRACE_LEVEL("EnzoSolverHg::restrict_recv",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverHg::restrict_recv()");
  
  if (enzo_block->mg_sync_next()) {
    //    LOCAL_SUM("3 B",ib_);
    begin_cycle_<T>(enzo_block);
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::solve_coarse(EnzoBlock * enzo_block) throw()
/// 
/// [*] coarse solve
///
///      solve A X = B
///      end_cycle()
{
  TRACE_LEVEL("EnzoSolverHg::solve_coarse",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverHg::solver_coarse()");
 
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

  /// Apply smoother for coarse-grid solve

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  Solver * solve_coarse = simulation->problem()->solver(index_solve_coarse_);

  solve_coarse->apply(A_,ix_,ib_,enzo_block);

  if (enzo_block->level() < max_level_) {

    prolong_send_<T>(enzo_block);

  } 
  end_cycle<T>(enzo_block);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::prolong_send_(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong send
///
///      for child           
///         pack X
///         child.prolong_recv(X)
{
  TRACE_LEVEL("EnzoSolverHg::prolong_send",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverHg::prolong_send()");

  ItChild it_child(enzo_block->rank());
  int ic3[3];
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

    // CkPrintf ("%s p_solver_hg_prolong_recv\n",
    // 	      enzo_block->name().c_str(), enzo_block->level());
    // fflush(stdout);
    CkCallback (CkIndex_EnzoBlock::p_solver_hg_prolong_recv<T>(NULL), 
		CkArrayIndexIndex(index_child),
		enzo_block->proxy_array()).send(field_message);

    delete [] array;

  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_hg_prolong_recv(FieldMsg * field_message)
/// [*]
{

  TRACE_MG (this,"EnzoBlock::p_solver_hg_prolong_recv()");
  
  // Unpack "C" vector data from children

  int if3[3] = {0,0,0};
  bool lg3[3] = {false,false,false};
  std::vector<int> field_list;
  field_list.push_back(data()->field().field_id("C"));

  // copy data from field_message to this EnzoBlock

  EnzoSolverHg * solver = 
    static_cast<EnzoSolverHg*> (this->solver());

  FieldFace * field_face = create_face 
    (if3, field_message->ic3, lg3, refresh_fine, field_list);

  field_face->set_prolong(solver->prolong());

  field_face->array_to_face (field_message->a, data()->field());

  delete field_face;

  // GET DATA
  
  //  delete field_message;
  delete field_message;

  solver -> prolong_recv<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::prolong_recv(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong recv
///
///      [ unpack C ]
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{

  TRACE_LEVEL("EnzoSolverHg::prolong_recv",enzo_block);
  TRACE_MG (enzo_block,"EnzoSolverHg::prolong_recv()");
  
  //  LOCAL_SUM("5 X",ix_);
  //  LOCAL_SUM("5 C",ic_);

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

  Refresh refresh (4,0,neighbor_level, sync_face,8);

  refresh.add_all_fields(enzo_block->data()->field().field_count());

  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_hg_post_smooth<T>(NULL),&refresh);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_hg_post_smooth(CkReductionMsg * msg)
/// [*]
{
  TRACE_MG (this,"EnzoBlock::p_solver_hg_post_smooth()");
  
  delete msg;

  EnzoSolverHg * solver = 
    static_cast<EnzoSolverHg*> (this->solver());

  solver -> post_smooth<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::post_smooth(EnzoBlock * enzo_block) throw()
/// post smooth
/// [*]
///
///      smooth.apply (A,X,B)
///      end_cycle()
{
  TRACE_LEVEL("EnzoSolverHg::post_smooth",enzo_block);
  TRACE_MG(enzo_block,"EnzoSolverHg::post_smooth()");

  //  smooth_post_->apply(A_,ix_,ib_,enzo_block);

  //  LOCAL_SUM("6 X",ix_);

  if (enzo_block->level() < max_level_) {

    prolong_send_<T>(enzo_block);

  } 

  end_cycle<T>(enzo_block);
}


//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_hg_end_cycle(CkReductionMsg * msg)
{
  /// EnzoBlock accumulates global contributions to SUM(B) and COUNT(B)
  EnzoSolverHg* solver = 
    static_cast<EnzoSolverHg*> (this->solver());
  
  long double* data = (long double*) msg->getData();

  // CkPrintf ("sum(X) = %Lf sum(C) = %Lf\n",data[0],data[1]);

#ifdef DEBUG_SOLVER_HG  
  CkPrintf ("%d %s DEBUG_SOLVER_HG bs %lf bc %lf\n",
	    __LINE__,name().c_str(),double(data[0]),double(data[1]));
#endif  

  delete msg;

  /// Start solver
  solver->end_cycle<T>(this);

}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverHg::end_cycle(EnzoBlock * enzo_block) throw()
/// 
/// [X] end cycle
///
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
{
  TRACE_LEVEL("EnzoSolverHg::end_cycle",enzo_block);
  Field field = enzo_block->data()->field();

  const int level = enzo_block->level();

  TRACE_MG(enzo_block,"EnzoSolverHg::end_cycle()");
  
  enzo_block->mg_iter_increment();

  if (is_converged_(enzo_block)) {

    //    exit_solver_<T>(enzo_block, return_converged);
    // implicitly exit solver

  }

  if (level >= max_level_) {

    begin_cycle_<T>(enzo_block);

  }
}


//======================================================================

bool EnzoSolverHg::is_converged_(EnzoBlock * enzo_block) const
/// [*]
{
  TRACE_MG(enzo_block,"EnzoSolverHg::is_converged");
  return (enzo_block->mg_iter() >= iter_max_);
}

//======================================================================
