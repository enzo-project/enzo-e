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
///      smoother.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
///
///  p_restrict_send(X)
///
///      A.residual(R,B,X)
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
///      smoother.apply (A,X,B)
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

#define MG_VERBOSE(FUNCTION)						\
  Monitor * monitor = enzo_block->simulation()->monitor();		\
  monitor->verbose(stdout,"Method", "%s "FUNCTION,enzo_block->name().c_str());

#define VERBOSE(FUNCTION)						\
  Monitor * monitor = simulation()->monitor();				\
  monitor->verbose(stdout,"Method", "%s "FUNCTION,name().c_str());

//----------------------------------------------------------------------

extern CkReduction::reducerType r_method_gravity_mg_type;

//----------------------------------------------------------------------

EnzoMethodGravityMg0::EnzoMethodGravityMg0 
(const FieldDescr * field_descr, 
 int rank,
 double grav_const, int iter_max, int monitor_iter,
 bool is_singular,
 Compute * smooth,
 Restrict * restrict,
 Prolong * prolong,
 int min_level,
 int max_level) 
  : Method(), 
    A_(new EnzoMatrixLaplace),
    smooth_(smooth),
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
    iter_(0),
    min_level_(min_level),
    max_level_(max_level)
  // [*]
{

  if (max_level_ > 0) {
    WARNING1 ("EnzoMethodGravityMg0::EnzoMethodGravityMg0()",
	     "Solver Mg0 is for root-level solves only: "
	     "changing max_level_ from %d to 0",max_level_);
    max_level_ = 0;
  }
  if (min_level_ >= 0) {
    ERROR1 ("EnzoMethodGravityMg0::EnzoMethodGravityMg0()",
	    "Solver Mg0 requires min_level = %d to be < 0",
	    min_level_);
  }
  /// Initialize default Refresh

  int ir = add_refresh(1,0,sync_barrier);
  refresh(ir)->add_all_fields(field_descr->field_count());

  int il = add_refresh(1,0,sync_face);
  refresh(il)->add_all_fields(field_descr->field_count());

  irho_ = field_descr->field_id("density");
  iphi_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  ir_ = field_descr->field_id("R");
  ix_ = field_descr->field_id("X");
  ic_ = field_descr->field_id("C");
}

//----------------------------------------------------------------------

EnzoMethodGravityMg0::~EnzoMethodGravityMg0 () throw()
// [*]
{
  delete smooth_;
  delete prolong_;
  delete restrict_;

  smooth_ = NULL;
  prolong_ = NULL;
  restrict_ = NULL;
}


//----------------------------------------------------------------------

void EnzoMethodGravityMg0::compute ( Block * block) throw()
// [*]
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  MG_VERBOSE("compute()");

  Field field = block->data()->field();
  precision_ = field.precision(irho_);

  if      (precision_ == precision_single)    
    enter_solver_<float>      (enzo_block);
  else if (precision_ == precision_double)    
    enter_solver_<double>     (enzo_block);
  else if (precision_ == precision_quadruple) 
    enter_solver_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityMg0()", "precision %d not recognized", precision_);
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
  
  iter_ = 0;

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

  int mx,my,mz;
  field.dimensions(irho_,&mx,&my,&mz);

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	int i = ix + mx*(iy + my*iz);
	B[i] = - 4.0 * (cello::pi) * grav_const_ * rho[i];
	X[i] = 0.0;
	R[i] = B[i];
	C[i] = 0.0;
      }
    }
  }

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

  if (is_converged_()) {

    exit_solver_<T>(enzo_block, return_converged);

  } else if (level == min_level_) {

    solve_coarse_<T>(enzo_block);

  } else {

    enzo_block->refresh_enter
      (CkIndex_EnzoBlock::p_mg0_pre_smooth<T>(NULL),refresh(1));

  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_pre_smooth(CkReductionMsg * msg)
/// [*]
{
  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // GET DATA
  
  //  delete msg;
  delete msg;

  method -> pre_smooth<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::pre_smooth(EnzoBlock * enzo_block) throw()
/// [*]
///
///      smoother.apply (A,X,B)
///      callback = p_restrict_send()
///      call refresh (X,level,"level")
{
  MG_VERBOSE("pre_smooth_()");

  smooth_->compute(enzo_block);

    enzo_block->refresh_enter
      (CkIndex_EnzoBlock::p_mg0_restrict_send<T>(NULL),refresh(1));
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_restrict_send(CkReductionMsg * msg)
/// [*]
{
  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // GET DATA
  
  //  delete msg;
  delete msg;

  method -> restrict_send<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::restrict_send(EnzoBlock * enzo_block) throw()
/// 
/// [ ] restrict send
///
///      A.residual(R,B,X)
///      index_parent.p_restrict_recv(R)
{

  MG_VERBOSE("restrict_send_()");

  A_->residual(ir_, ib_, ix_, enzo_block);

  Index index_parent = enzo_block->index().index_parent(min_level_);
  
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_restrict_recv(CkReductionMsg * msg)
/// [*]
{
  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // GET DATA
  
  //  delete msg;
  delete msg;

  method -> restrict_recv<T>(this);
}


//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::restrict_recv(EnzoBlock * enzo_block) throw()
/// 
/// [ ] restrict recv
///
///      unpack B
///      if (sync.next())
///          begin_cycle()
{
  MG_VERBOSE("restrict_recv_()");
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::solve_coarse_(EnzoBlock * enzo_block) throw()
/// 
/// [ ] coarse solve
///
///      solve A X = B
///      end_cycle()
{
  MG_VERBOSE("solve_coarse_()");
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::prolong_send_(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong send
///
///      for child           
///         child.prolong_recv(X)
{
  MG_VERBOSE("prolong_send_()");
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_prolong_recv(CkReductionMsg * msg)
/// [*]
{
  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // GET DATA
  
  //  delete msg;
  delete msg;

  method -> prolong_recv<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::prolong_recv(EnzoBlock * enzo_block) throw()
/// 
/// [ ] prolong recv
///
///      unpack C         
///      X = X + C
///      callback = p_post_smooth()
///      call refresh (X,"level")
{
  MG_VERBOSE("prolong_recv_()");
}

//----------------------------------------------------------------------

template <class T>
void EnzoBlock::p_mg0_post_smooth(CkReductionMsg * msg)
/// [*]
{
  EnzoMethodGravityMg0 * method = 
    static_cast<EnzoMethodGravityMg0*> (this->method());

  // GET DATA
  
  //  delete msg;
  delete msg;

  method -> post_smooth<T>(this);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::post_smooth(EnzoBlock * enzo_block) throw()
/// post smooth
/// [*]
///
///      smoother.apply (A,X,B)
///      end_cycle()
{
  MG_VERBOSE("post_smooth_()");

  smooth_->compute(enzo_block);

  end_cycle_<T>(enzo_block);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg0::end_cycle_(EnzoBlock * enzo_block) throw()
/// 
/// [ ] end cycle
///
///      ++iter
///      if (level < max_level)
///         prolong_send(X)
///      else
///         begin_cycle()
{
  MG_VERBOSE("end_cycle_()");
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

  int mx,my,mz;
  field.dimensions(irho_,&mx,&my,&mz);

  copy_(phi,X,mx,my,mz,enzo_block->is_leaf());

  FieldDescr * field_descr = field.field_descr();

  EnzoComputeAcceleration compute_acceleration (field_descr,rank_, true,2);

  compute_acceleration.compute(enzo_block);

  monitor_output_ (enzo_block);

  enzo_block->compute_done();
}

//======================================================================

void EnzoMethodGravityMg0::monitor_output_(EnzoBlock * enzo_block) throw()
/// [*]
{
  if (enzo_block->index().is_root()) {

    Monitor * monitor = enzo_block->simulation()->monitor();

    monitor->print ("Enzo", "Mg0 iter %04d  rr %g",iter_, (double)(rr_/rr0_));
  }
}

//======================================================================

bool EnzoMethodGravityMg0::is_converged_() const
/// [*]
{
  return (iter_ >= iter_max_);
}

//======================================================================
