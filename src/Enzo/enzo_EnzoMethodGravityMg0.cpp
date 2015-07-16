// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityMg0.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityMg0 class
///
/// Multigrid method on a non-adaptive mesh.
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
/// - phi                        computed gravitational potential
/// - rho                        density field
/// - acceleration_x             acceleration along X-axis
/// - acceleration_y (rank >= 2) acceleration along Y-axis
/// - acceleration_z (rank >= 3) acceleration along Z-axis

#include "cello.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#define DEBUG_MG
/* #define DEBUG_REFRESH_ORDER */



#ifdef DEBUG_MG
#define MG_VERBOSE(FUNCTION)						\
  Monitor * monitor = enzo_block->simulation()->monitor();		\
  monitor->verbose(stdout,"Method", "%s "FUNCTION,enzo_block->name().c_str());

#define VERBOSE(FUNCTION)						\
  Monitor * monitor = simulation()->monitor();				\
  monitor->verbose(stdout,"Method", "%s "FUNCTION,name().c_str());
#define DEBUG_X							\
  {								\
    Data * data = enzo_block->data();				\
    Field field = data->field();				\
								\
    T * X = (T*) field.values(ix_);				\
    T * B = (T*) field.values(ib_);				\
    std::string name = enzo_block->name();			\
    if (name == "B0_0") {					\
      CkPrintf ("%s:%d DEBUG %s X(1,[19:21]) = [%g %g %g]\n",	\
	      __FILE__,__LINE__,name.c_str(),			\
	      X[43],X[44],X[45]);				\
      CkPrintf ("%s:%d DEBUG %s B(1,[19:21]) = [%g %g %g]\n",	\
	      __FILE__,__LINE__,name.c_str(),			\
	      B[43],B[44],B[45]);				\
    } else if (name == "B1_0") {				\
      CkPrintf ("%s:%d DEBUG %s X(1,[19:21]) = [%g %g %g]\n",	\
	      __FILE__,__LINE__,name.c_str(),			\
	      X[43],X[44],X[45]);				\
      CkPrintf ("%s:%d DEBUG %s B(1,[19:21]) = [%g %g %g]\n",	\
	      __FILE__,__LINE__,name.c_str(),			\
	      B[43],B[44],B[45]);				\
    }								\
  }
#else
#   define DEBUG_X	;
#   define MG_VERBOSE(X) /* */ ;
#   define VERBOSE(X) /* */ ;
#endif

#ifdef DEBUG_MG
#   define PRINT_BDRX				\
  for (int ix=0; ix<mx; ix++) {			\
    for (int iy=0; iy<my; iy++) {		\
      for (int iz=0; iz<mz; iz++) {		\
	int i = ix + mx*(iy + my*iz);		\
	CkPrintf ("[%d %d] %d B=%6.3g D=%6.3g R=%6.3g X=%6.3g\n",	\
		ix,iy,i,B[i],D[i],R[i],X[i]);	\
      }						\
    }						\
  }
#else 
#  define PRINT_BDRX ;
#endif

#define PRINT_X \
    printf ("%s:%d DEBUG_X %s X[4,0] = %g  X[4,20]\n", \
	    __FILE__,__LINE__,block->name().c_str(),X[4]);

//----------------------------------------------------------------------

extern CkReduction::reducerType r_method_gravity_mg_type;

//----------------------------------------------------------------------

EnzoMethodGravityMg0::EnzoMethodGravityMg0 
(const FieldDescr * field_descr, 
 int rank,
 double grav_const, int iter_max, int monitor_iter, 
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
  : Method(), 
    A_(new EnzoMatrixLaplace),
    smooth_pre_(NULL),
    smooth_coarse_(NULL),
    smooth_post_(NULL),
    restrict_(restrict),
    prolong_(prolong),
    is_singular_(is_singular),
    rank_(rank),
    grav_const_(grav_const),
    iter_max_(iter_max), 
    monitor_iter_(monitor_iter),
    rr_(0),rr0_(0),
    irho_(0),  iphi_(0),
    ib_(0), ir_(0), ix_(0), ic_(0),
    min_level_(min_level),
    max_level_(max_level),
    mx_(0),my_(0),mz_(0)

  // [*]
{

  if (max_level_ > 0) {
    WARNING1 ("EnzoMethodGravityMg0::EnzoMethodGravityMg0()",
	     "Solver Mg0 is for root-level solves only: "
	     "changing max_level_ from %d to 0",max_level_);
    max_level_ = 0;
  }
  if (min_level_ > 0) {
    ERROR1 ("EnzoMethodGravityMg0::EnzoMethodGravityMg0()",
	    "Solver Mg0 requires min_level = %d to be <= 0",
	    min_level_);
  }

  /// Initialize default Refresh

  add_refresh(4,0,neighbor_level,sync_barrier);
  refresh(0)->add_all_fields(field_descr->field_count());

  irho_ = field_descr->field_id("density");
  iphi_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  ir_ = field_descr->field_id("R");
  ix_ = field_descr->field_id("X");
  ic_ = field_descr->field_id("C");
  id_ = field_descr->field_id("D");
  if (smooth == "jacobi") {
    smooth_pre_ = new EnzoComputeSmoothJacobi 
      (A_,ix_,ib_,ir_,id_, smooth_weight, smooth_count_pre);
    smooth_coarse_ = new EnzoComputeSmoothJacobi 
      (A_,ix_,ib_,ir_,id_, smooth_weight, smooth_count_coarse);
    smooth_post_ = new EnzoComputeSmoothJacobi 
      (A_,ix_,ib_,ir_,id_, smooth_weight, smooth_count_post);

  }
}

//----------------------------------------------------------------------

EnzoMethodGravityMg0::~EnzoMethodGravityMg0 () throw()
// [*]
{
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

void EnzoMethodGravityMg0::compute ( Block * block) throw()
// [*]
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  MG_VERBOSE("compute()");

  const int level = block->level();
  if (level < min_level_ || level > max_level_) {
    block->compute_done();
  }
  Field field = enzo_block->data()->field();

  field.dimensions(irho_,&mx_,&my_,&mz_);

  // Initialize child counter for restrict synchronization
  enzo_block->mg_sync_set_stop(NUM_CHILDREN(enzo_block->rank()));
  enzo_block->mg_sync_reset();

  int precision = field.precision(irho_);

  if      (precision == precision_single)    
    enter_solver_<float>      (enzo_block);
  else if (precision == precision_double)    
    enter_solver_<double>     (enzo_block);
  else if (precision == precision_quadruple) 
    enter_solver_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityMg0()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::enter_solver_ (EnzoBlock * enzo_block) throw()
/// [*]
///
///     iter = 0
///     initialize X,B,R,C
///     if (level == max_level) 
///        begin_cycle()
{

  MG_VERBOSE("enter_solver_()");
  
  enzo_block->mg_iter_clear();

  const int level = enzo_block->level();

  Data * data = enzo_block->data();
  Field field = data->field();

  T * rho = (T*) field.values(irho_);
    
  T * B = (T*) field.values(ib_);
  T * X = (T*) field.values(ix_);
  T * R = (T*) field.values(ir_);
  T * C = (T*) field.values(ic_);

  // X = 0
  // B = -h^2 * 4 * PI * G * rho
  // R = B ( residual with X = 0 )
  // C = 0

  int gx,gy,gz;
  int mx,my,mz;
  field.dimensions  (0,&mx,&my,&mz);
  field.ghost_depth (0,&gx,&gy,&gz);

  const int ix0 = 0;
  const int iy0 = 0;
  const int iz0 = 0;
  for (int iz=iz0; iz<mz-iz0; iz++) {
    for (int iy=iy0; iy<my-iy0; iy++) {
      for (int ix=ix0; ix<mx-ix0; ix++) {
	int i = ix + mx*(iy + my*iz);
	B[i] = - 4.0 * (cello::pi) * grav_const_ * rho[i];
	X[i] = 0.0;
	R[i] = B[i];
	C[i] = 0.0;
      }
    }
  }

  DEBUG_X;
  // start the MG V-cycle with the root level
  if (enzo_block->level() == max_level_)
    begin_cycle_<T>(enzo_block);

}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::begin_cycle_(EnzoBlock * enzo_block) throw()
/// [*]
///
///     if (converged()) exit()
///     if (level == min_level) then
///        coarse_solve(A,X,B)
///     else
///        callback = p_pre_smooth()
///        call refresh (X,"level")
{
  MG_VERBOSE("begin_cycle_()");

  const int level = enzo_block->level();

  if (is_converged_(enzo_block)) {

    exit_solver_<T>(enzo_block, return_converged);

  } else if (level == min_level_) {
    // TODO REFRESH X
    solve_coarse_<T>(enzo_block);

  } else {

    const Data * data = enzo_block->data();
    const FieldDescr * field_descr = data->field_descr();
    const int num_fields = field_descr->field_count();

    Refresh refresh (4,0,neighbor_level, sync_face);
    // refresh.add_all_fields (enzo_block->data()->field_descr()->field_count());
    refresh.add_field (ix_);

    DEBUG_X;
    
    enzo_block->refresh_enter
      (CkIndex_EnzoBlock::p_mg0_pre_smooth<T>(NULL),&refresh);
    // // Skip refresh
    // pre_smooth<T>(enzo_block);

  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_pre_smooth(CkReductionMsg * msg)
/// [*]
{
  VERBOSE("p_mg0_pre_smooth()");
  delete msg;

  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  method -> pre_smooth<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::pre_smooth(EnzoBlock * enzo_block) throw()
/// [*]
///
///      smooth.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  MG_VERBOSE("pre_smooth_()");

  Data * data = enzo_block->data();
  Field field = data->field();
  T * X = (T*) field.values(ix_);
  DEBUG_X;
  smooth_pre_->compute(enzo_block);
  DEBUG_X;

  Refresh refresh (4,0,neighbor_level, sync_face);
  //  refresh.add_all_fields (enzo_block->data()->field_descr()->field_count());
    refresh.add_field (ix_);

  refresh.set_active(true);
  refresh.print();
  DEBUG_X;

  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_mg0_restrict_send<T>(NULL),&refresh);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // // Skip refresh
  // restrict_send<T>(enzo_block);


}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_restrict_send(CkReductionMsg * msg)
/// [*]
{
  VERBOSE("p_mg0_restrict_send()");
  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // GET DATA
  
  delete msg;

  method -> restrict_send<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::restrict_send(EnzoBlock * enzo_block) throw()
/// 
/// [*] restrict send
///
///      A.residual(R,B,X)
///      pack R
///      index_parent.p_restrict_recv(R)
{
  Data * data = enzo_block->data();
  Field field = data->field();
  T * X = (T*) field.values(ix_);

  DEBUG_X;
  MG_VERBOSE("restrict_send_()");

  A_->residual(ir_, ib_, ix_, enzo_block);

  DEBUG_X;

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

  int narray; 
  char * array;
  FieldFace * field_face = 
    enzo_block->load_face(&narray,&array, if3, ic3, lg3, 
			  op_array_restrict, field_list);

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
  CkCallback (CkIndex_EnzoBlock::p_mg0_restrict_recv<T>(NULL), 
	      CkArrayIndexIndex(index_parent),
	      enzo_block->proxy_array()).send(field_message);

  delete field_face;

}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_restrict_recv(FieldMsg * field_message)
/// [*]
{

  VERBOSE("p_mg0_restrict_recv()");

  // Unpack "B" vector data from children

  // Face (0,0,0) is the entire block
  int if3[3] = {0,0,0};
  // exclude ghost zones
  bool lg3[3] = {false,false,false};
  // receive in vector "B"
  std::vector<int> field_list;
  field_list.push_back(data()->field().field_id("B"));

  // copy data from field_message to this EnzoBlock

  store_face_(field_message->n,
	      field_message->a,   if3, 
	      field_message->ic3, lg3,
	      op_array_restrict,
	      field_list);

  delete field_message;

  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // continue with EnzoMethodGravityMg0
  method -> restrict_recv<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::restrict_recv(EnzoBlock * enzo_block) throw()
/// 
/// [*] restrict recv
///
///      [ unpack B ]
///      if (sync.next())
///          begin_cycle()
{
  MG_VERBOSE("restrict_recv_()");

  if (enzo_block->mg_sync_next()) {
    begin_cycle_<T>(enzo_block);
  }
  
  
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::solve_coarse_(EnzoBlock * enzo_block) throw()
/// 
/// [*] coarse solve
///
///      solve A X = B
///      end_cycle()
{
  MG_VERBOSE("solve_coarse_()");

  /// Apply smoother
  smooth_coarse_->compute(enzo_block);

  if (enzo_block->level() < max_level_) {

    prolong_send_<T>(enzo_block);

  } 
  end_cycle_<T>(enzo_block);

}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::prolong_send_(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong send
///
///      for child           
///         pack X
///         child.prolong_recv(X)
{
  MG_VERBOSE("prolong_send_()");
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

    int narray; 
    char * array;
    FieldFace * field_face = 
      enzo_block->load_face(&narray,&array, if3, ic3, lg3, 
			    op_array_prolong, field_list);

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

    CkCallback (CkIndex_EnzoBlock::p_mg0_prolong_recv<T>(NULL), 
		CkArrayIndexIndex(index_child),
		enzo_block->proxy_array()).send(field_message);
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_prolong_recv(FieldMsg * field_message)
/// [*]
{
  VERBOSE("p_mg0_prolong_recv()");

  // Unpack "C" vector data from children

  // Face (0,0,0) is the entire block
  int if3[3] = {0,0,0};
  // exclude ghost zones
  bool lg3[3] = {false,false,false};
  // receive in vector "C"
  std::vector<int> field_list;
  field_list.push_back(data()->field().field_id("C"));

  // copy data from field_message to this EnzoBlock

  store_face_(field_message->n,
	      field_message->a,   if3, 
	      field_message->ic3, lg3,
	      op_array_prolong,
	      field_list);

  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // GET DATA
  
  //  delete field_message;
  delete field_message;

  method -> prolong_recv<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::prolong_recv(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong recv
///
///      [ unpack C ]
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{
  MG_VERBOSE("prolong_recv_()");

  DEBUG_X;

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X = (T*) field.values(ix_);
  T * C = (T*) field.values(ic_);

  zaxpy_(X,1.0,X,C);

  DEBUG_X;

  Refresh refresh (4,0,neighbor_level, sync_face);
  // refresh.add_all_fields (enzo_block->data()->field_descr()->field_count());
    refresh.add_field (ix_);

  enzo_block->refresh_enter
    (CkIndex_EnzoBlock::p_mg0_post_smooth<T>(NULL),&refresh);
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_post_smooth(CkReductionMsg * msg)
/// [*]
{
  VERBOSE("p_mg0_post_smooth()");

  //  delete msg;
  delete msg;

  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  method -> post_smooth<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::post_smooth(EnzoBlock * enzo_block) throw()
/// post smooth
/// [*]
///
///      smooth.apply (A,X,B)
///      end_cycle()
{
  MG_VERBOSE("post_smooth_()");

  smooth_post_->compute(enzo_block);

  if (enzo_block->level() < max_level_) {

    prolong_send_<T>(enzo_block);

  } 

  end_cycle_<T>(enzo_block);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::end_cycle_(EnzoBlock * enzo_block) throw()
/// 
/// [X] end cycle
///
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
{
  MG_VERBOSE("end_cycle_()");

  enzo_block->mg_iter_increment();

  if (is_converged_(enzo_block)) {

    exit_solver_<T>(enzo_block, return_converged);

  }

  if (enzo_block->level() >= max_level_) {

    begin_cycle_<T>(enzo_block);

  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::exit_solver_ 
( EnzoBlock * enzo_block, int retval ) throw ()
/// 
/// [ ] exit solver
///
///      acceleration.compute(X)
{
  MG_VERBOSE("exit_solver_()");

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X         = (T*) field.values(ix_);
  T * phi = (T*) field.values(iphi_);

  DEBUG_X;

  int mx,my,mz;
  field.dimensions(irho_,&mx,&my,&mz);

  copy_(phi,X,mx,my,mz,enzo_block->is_leaf());

  DEBUG_X;

  FieldDescr * field_descr = field.field_descr();

  EnzoComputeAcceleration compute_acceleration (field_descr,rank_, true,2);

  compute_acceleration.compute(enzo_block);

  monitor_output_ (enzo_block);

  DEBUG_X;
  enzo_block->compute_done();
}

//======================================================================

void EnzoMethodGravityMg0::monitor_output_(EnzoBlock * enzo_block) throw()
/// [*]
{
  if (enzo_block->index().is_root()) {

    Monitor * monitor = enzo_block->simulation()->monitor();

    monitor->print ("Enzo", "Mg0 iter %04d  rr %g",
		    enzo_block->mg_iter(),
		    (double)(rr_/rr0_));
  }
}

//======================================================================

bool EnzoMethodGravityMg0::is_converged_(EnzoBlock * enzo_block) const
/// [*]
{
  return (enzo_block->mg_iter() >= iter_max_);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::zaxpy_ 
(T * Z, double a, const T * X, const T * Y) const throw()
{
  for (int iz=0; iz<mz_; iz++) {
    for (int iy=0; iy<my_; iy++) {
      for (int ix=0; ix<mx_; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	Z[i] = a * X[i] + Y[i];
      }
    }
  }
}

//======================================================================
