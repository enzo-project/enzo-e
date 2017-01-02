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
///    if (level == max_level)
///       if (converged()) exit()
///    if (level == min_level) then
///       coarse_solve()     solve $A_h X_h = B_h$
///    else
///       p_pre_smooth()     smooth $A_h X_h = B_h$
///       p_residual()       $R_h = B_h - A_h * X_h$
///       p_restrict ()      $B_H = I_h^H R_h$
///       MG()               solve $A_H X_H = B_H$  (repeat for W-cycle)
///       p_prolong ()       $X_h = X_h + I_H^h X_H$
///       p_post_smooth()    smooth $A_h X_h = B_h$
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
///     initialize X,B,R,C
///     if (level == max_level) 
///        begin_cycle()
///
///  begin_cycle()
///
///     if (converged()) exit()
///     if (level == min_level) then
///        coarse_solve(A,X,B)
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
///      A.residual(R,B,X)
///      pack R
///      index_parent.p_restrict_recv(R)
///
///  p_restrict_recv(B)
///
///      unpack B
///      if (sync.next())
///          begin_cycle()
///      
///  coarse_solve(A,X,B)
///
///      solve A X = B
///      end_cycle()
///
///  prolong_send(X)
///
///      for child           
///         pack X
///         child.prolong_recv(X)
///
///  prolong_recv(C)
///
///      unpack C         
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
///
///  p_post_smooth(A,X,B)
/// 
///      smooth.apply (A,X,B)
///      end_cycle()
///
///  end_cycle()
///
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
///
///  exit_solver()
///  
///      acceleration.compute(X)
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
/// - acceleration_x             acceleration along X-axis
/// - acceleration_y (rank >= 2) acceleration along Y-axis
/// - acceleration_z (rank >= 3) acceleration along Z-axis

#include "cello.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

// #define DEBUG_MG

#ifdef DEBUG_MG
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

// //----------------------------------------------------------------------

// extern CkReduction::reducerType r_method_gravity_mg_type;

//----------------------------------------------------------------------

EnzoSolverMg0::EnzoSolverMg0 
(const FieldDescr * field_descr, 
 int monitor_iter,
 int rank,
 int iter_max,
 std::string smooth,
 double      smooth_weight,
 int         smooth_count_pre,
 int         smooth_count_coarse,
 int         smooth_count_post,
 bool is_singular,
 Restrict * restrict,
 Prolong * prolong,
 int min_level,
 int max_level) 
  : Solver(), 
    A_(new EnzoMatrixLaplace),
    smooth_pre_(NULL),
    smooth_coarse_(NULL),
    smooth_post_(NULL),
    restrict_(restrict),
    prolong_(prolong),
    is_singular_(is_singular),
    rank_(rank),
    iter_max_(iter_max), 
    monitor_iter_(monitor_iter),
    ib_(0), ic_(0), ir_(0), ix_(0),
    min_level_(min_level),
    max_level_(max_level),
    mx_(0),my_(0),mz_(0),
    nx_(0),ny_(0),nz_(0),
    gx_(0),gy_(0),gz_(0)

    // [*]
{
  EnzoBlock * null_block = NULL;
  TRACE_MG(null_block,"EnzoSolverMg0()");
  
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
  id_ = field_descr->field_id("D");

  if (smooth == "jacobi") {
    smooth_pre_ = new EnzoSolverJacobi
      (id_, ir_, smooth_weight, smooth_count_pre);
    smooth_coarse_ = new EnzoSolverJacobi
      (id_, ir_, smooth_weight, smooth_count_coarse);
    smooth_post_ = new EnzoSolverJacobi
      (id_, ir_, smooth_weight, smooth_count_post);
  }

}

//----------------------------------------------------------------------

EnzoSolverMg0::~EnzoSolverMg0 () throw()
// [*]
{
  EnzoBlock * null_block = NULL;
  TRACE_MG(null_block,"~EnzoSolverMg0()");

  delete smooth_pre_;
  delete smooth_coarse_;
  delete smooth_post_;
  delete prolong_;
  delete restrict_;

  smooth_pre_ = NULL;
  smooth_coarse_ = NULL;
  smooth_post_ = NULL;
  prolong_ = NULL;
  restrict_ = NULL;
}


//----------------------------------------------------------------------

void EnzoSolverMg0::apply
( Matrix * A, int ix, int ib, Block * block) throw()
{
  TRACE_MG(block,"EnzoSolverMg0::apply()");

  A_ = A;
  ix_ = ix;
  ib_ = ib;

  block->set_solver(this);
  Field field = block->data()->field();

  field.size           (&nx_,&ny_,&nz_);
  field.dimensions (ib_,&mx_,&my_,&mz_);
  field.ghost_depth(ib_,&gx_,&gy_,&gz_);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  const int level = block->level();

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
  TRACE_MG(enzo_block,"EnzoSolverMg0::enter_solver()");

  enzo_block->mg_iter_clear();

  Data * data = enzo_block->data();
  Field field = data->field();

  T * B = (T*) field.values(ib_);
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
	R[i] = B[i];
	C[i] = 0.0;
      }
    }
  }

  // start the MG V-cycle with the root level
  if (enzo_block->level() == max_level_) {
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
  TRACE_MG(enzo_block,"EnzoSolverMg0::begin_cycle()");
  const int level = enzo_block->level();

  if (is_converged_(enzo_block)) {

    TRACE_MG(enzo_block,"converged");
    // implicitly exit solver

  } else if (level == min_level_) {
    // TODO REFRESH X
    TRACE_MG(enzo_block,"calling coarse solve");
    solve_coarse_<T>(enzo_block);

  } else {

    TRACE_MG(enzo_block,"calling smoother");

    const Data * data = enzo_block->data();

    Refresh refresh (4,0,neighbor_level, sync_face);
    refresh.add_all_fields(enzo_block->data()->field().field_count());

    refresh.set_active(true);

    //    enzo_block->set_solver(this);
    enzo_block->refresh_enter
      (CkIndex_EnzoBlock::p_solver_mg0_pre_smooth(),&refresh);
    // // Skip refresh
    // pre_smooth<T>(enzo_block);

  }
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
void EnzoSolverMg0::pre_smooth(EnzoBlock * enzo_block) throw()
/// [*]
///
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  TRACE_MG(enzo_block,"EnzoSolverMg0::pre_smooth()");

  Data * data = enzo_block->data();
  Field field = data->field();
  T * X = (T*) field.values(ix_);

  smooth_pre_->apply(A_,ix_,ib_,enzo_block);

  Refresh refresh (4,0,neighbor_level, sync_face);
  refresh.add_all_fields(enzo_block->data()->field().field_count());

  refresh.set_active(true);

  //  enzo_block->set_solver(this);
  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_mg0_restrict_send<T>(NULL),&refresh);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // // Skip refresh
  // restrict_send<T>(enzo_block);


}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_mg0_restrict_send(CkReductionMsg * msg)
/// [*]
{
  TRACE_MG(this,"EnzoBlock::p_solver_mg0_restrict_send()");

  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  // GET DATA
  
  delete msg;

  solver -> restrict_send<T>(this);
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
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_send()");

  Data * data = enzo_block->data();
  Field field = data->field();

  A_->residual(ir_, ib_, ix_, enzo_block);


  Index index        = enzo_block->index();
  Index index_parent = index.index_parent(min_level_);
  int level          = index.level();
  
  // copy face data to FieldFace

  // face (0,0,0) is the entire block
  int if3[3] = {0,0,0};
  // find which child this block is in its parent
  int ic3[3];
  index.child(level,&ic3[0],&ic3[1],&ic3[2],min_level_);
  // exclude ghost zones
  bool lg3[3] = {false,false,false};
  // send vector "R" data only
  std::vector<int> field_list;
  field_list.push_back(ir_);

  // copy data from EnzoBlock to array via FieldFace

  FieldFace * field_face = 
    enzo_block->create_face (if3, ic3, lg3, 
			     refresh_coarse, field_list);
  int narray; 
  char * array;

  field_face->face_to_array(enzo_block->data()->field(),&narray,&array);

  delete field_face;

  // Create a FieldMsg for sending data to parent
  // (note: charm messages not deleted on send; are deleted on receive)

  FieldMsg * field_message  = new (narray) FieldMsg;
 
  /// WARNING: double copy

  // Copy FieldFace data to field_message

  // array size
  field_message->n = narray;
  // array values
  memcpy (field_message->a, array, narray);
  // child index
  field_message->ic3[0] = ic3[0];
  field_message->ic3[1] = ic3[1];
  field_message->ic3[2] = ic3[2];

  //  </COMMON CODE>

  // Send field_message to parent
  //  enzo_block->set_solver(this);

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

  // Unpack "B" vector data from children

  // Face (0,0,0) is the entire block
  int if3[3] = {0,0,0};
  // exclude ghost zones
  bool lg3[3] = {false,false,false};
  // receive in vector "B"
  std::vector<int> field_list;
  field_list.push_back(data()->field().field_id("B"));

  // copy data from field_message to this EnzoBlock


  int * ic3 = field_message->ic3;

  FieldFace * field_face = create_face 
    (if3, ic3, lg3, refresh_coarse, field_list);

  char * a = field_message->a;
  field_face->array_to_face(a, data()->field());
  delete field_face;

  delete field_message;

  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  // continue with EnzoSolverMg0

  solver -> restrict_recv<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::restrict_recv(EnzoBlock * enzo_block) throw()
/// 
/// [*] restrict recv
///
///      [ unpack B ]
///      if (sync.next())
///          begin_cycle()
{
  TRACE_MG(enzo_block,"EnzoSolverMg0::restrict_recv()");
  
  if (enzo_block->mg_sync_next()) {
    begin_cycle_<T>(enzo_block);
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::solve_coarse_(EnzoBlock * enzo_block) throw()
/// 
/// [*] coarse solve
///
///      solve A X = B
///      end_cycle()
{
  TRACE_MG(enzo_block,"EnzoSolverMg0::solver_coarse()");

  /// Apply smoother
  smooth_coarse_->apply(A_,ix_,ib_,enzo_block);

  if (enzo_block->level() < max_level_) {

    prolong_send_<T>(enzo_block);

  } 
  end_cycle_<T>(enzo_block);
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
  TRACE_MG(enzo_block,"EnzoSolverMg0::prolong_send()");

  ItChild it_child(enzo_block->rank());
  int ic3[3];
  while (it_child.next(ic3)) {

    Index index_child = enzo_block->index().index_child(ic3,min_level_);

    // face (0,0,0) is the entire block
    int if3[3] = {0,0,0};
    // exclude ghost zones
    bool lg3[3] = {false,false,false};
    // send vector "X" data only
    std::vector<int> field_list;
    field_list.push_back(ix_);

    // copy data from EnzoBlock to array via FieldFace

    // <COMMON CODE> in restrict_send_() and prolong_send_()

    FieldFace * field_face = enzo_block->create_face
      (if3, ic3, lg3, refresh_fine, field_list);

    int narray; 
    char * array;
    
    field_face->face_to_array (enzo_block->data()->field(),&narray,&array);

    delete field_face;

    // Create a FieldMsg for sending data to parent
    // (note: charm messages not deleted on send; are deleted on receive)
    FieldMsg * field_message  = new (narray) FieldMsg;

    /// WARNING: double copy

    // Copy FieldFace data to field_message

    // array size
    field_message->n = narray;
    // array values
    memcpy (field_message->a, array, narray);
    // child index
    field_message->ic3[0] = ic3[0];
    field_message->ic3[1] = ic3[1];
    field_message->ic3[2] = ic3[2];

    //  </COMMON CODE>

    CkCallback (CkIndex_EnzoBlock::p_solver_mg0_prolong_recv<T>(NULL), 
		CkArrayIndexIndex(index_child),
		enzo_block->proxy_array()).send(field_message);
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg * field_message)
/// [*]
{

  TRACE_MG (this,"EnzoBlock::p_solver_mg0_prolong_recv()");
  
  // Unpack "C" vector data from children

  // Face (0,0,0) is the entire block
  int if3[3] = {0,0,0};
  // exclude ghost zones
  bool lg3[3] = {false,false,false};
  // receive in vector "C"
  std::vector<int> field_list;
  field_list.push_back(data()->field().field_id("C"));

  // copy data from field_message to this EnzoBlock

  int n = field_message->n;
  char * a = field_message->a;
  int * ic3 = field_message->ic3;

  FieldFace * field_face = create_face 
    (if3, field_message->ic3, lg3, refresh_fine, field_list);

  field_face->array_to_face (field_message->a, data()->field());

  delete field_face;

  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  // GET DATA
  
  //  delete field_message;
  delete field_message;

  solver -> prolong_recv<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::prolong_recv(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong recv
///
///      [ unpack C ]
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{

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

  Refresh refresh (4,0,neighbor_level, sync_face);

  refresh.add_all_fields(enzo_block->data()->field().field_count());

  //  enzo_block->set_solver(this);
  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_solver_mg0_post_smooth<T>(NULL),&refresh);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_solver_mg0_post_smooth(CkReductionMsg * msg)
/// [*]
{
  TRACE_MG (this,"EnzoBlock::p_solver_mg0_post_smooth()");
  
  delete msg;

  EnzoSolverMg0 * solver = 
    static_cast<EnzoSolverMg0*> (this->solver());

  solver -> post_smooth<T>(this);
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
  TRACE_MG(enzo_block,"EnzoSolverMg0::post_smooth()");

  smooth_post_->apply(A_,ix_,ib_,enzo_block);

  if (enzo_block->level() < max_level_) {

    prolong_send_<T>(enzo_block);

  } 

  end_cycle_<T>(enzo_block);
}

//----------------------------------------------------------------------

template <class T>
void EnzoSolverMg0::end_cycle_(EnzoBlock * enzo_block) throw()
/// 
/// [X] end cycle
///
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
{
  TRACE_MG(enzo_block,"EnzoSolverMg0::end_cycle()");
  
  enzo_block->mg_iter_increment();

  if (is_converged_(enzo_block)) {

    //    exit_solver_<T>(enzo_block, return_converged);
    // implicitly exit solver

  }

  if (enzo_block->level() >= max_level_) {

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
