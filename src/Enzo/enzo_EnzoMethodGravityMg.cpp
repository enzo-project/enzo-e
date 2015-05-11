// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityMg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityMg class
///
/// Multigrid method on an adaptive mesh.  We use the MLAT (Multilevel
/// Adaptive Technique) of Brandt '77, based on the FAS formulation of
/// Multigrid.
///
///======================================================================
///
///  "Coarse" view of Multigrid solver
///
///  MG(A,X,B)
///
///    enter_solver()         initialize solver
///    begin_cycle()          initialize a V-cycle
///    send_faces()           send face data to neighbors
///  p_receive_face()            receive face daat from neighbors
///    pre_smooth()           apply the smoother
///  p_restrict()             receive restricted data from parent
///    compute_b()            compute the right-hand side B
///    coarse_solve()         solve the coarse-grid equation
///    p_end_coarse_solve()   callback after the coarse solver
///  p_prolong()              receive prolonged data from children
///    update_x()             update the solution X
///    post_smooth()          apply the smoother
///    end_cycle()            finalize V-cycle
///    exit_solver()          finalize solver
///
///
/// ======================================================================
///
///  "Fine" view of Multigrid algorithm
///
///  --------------------
///  enter_solver()
///
///     initialize iter_
///     call begin_cycle()
///
///  --------------------
///  begin_cycle()
///
///    if converged, then return
///    call send_faces()
///
///  --------------------
///  send_faces()
///
///    if leaf block, then
///       for each neighbor
///           pack face data
///           remote call p_receive_face() on neighbor
///       for each face
///           if neighbor in same level exists,
///                and it is not a leaf
///              pack face data
///              remote call p_receive_face() on neighbor 
///
///  --------------------
///  p_receive_face()
///
///     unpack ghost data
///     if all ghost data available
///        if coarsest level
///            call coarse_solve()
///        else 
///            if not leaf block and received restrict
///                call compute_b()
///            else
///                call pre_smooth()
///
///  --------------------
///  pre_smooth()
///
///      apply the smoother to A X = B
///      Y = A*X [ need refreshed X? ]
///      pack B, X, Y fields
///      remote call p_restrict() on parent
///
///  --------------------
///  p_restrict(B,X,Y)
///
///      unpack vectors B,X,Y
///      call send_faces()
///      if all ghost data available
///          if coarsest level
///             call coarse_solve()
///          else
///             call compute_b()
///
///  --------------------
///  compute_b()
///
///     Compute right hand side B:
///         B = B - Y
///         Y = A*X [ need refreshed X? ]
///         B = B + Y
///     call pre_smooth()
///
///  --------------------
///  begin_coarse_solve()
///     
///      [ need refreshed X?]
///      initiate the coarse solver with p_end_coarse_solve() callback
///
///  --------------------
///  p_end_coarse_solve()
///
///      if not leaf node
///          for each child
///              pack correction X
///              remote call p_prolong() on child
///      else
///          call update_x()
///
///  --------------------
///  p_prolong()
///
///     unpack correction Y
///     call update_x()
///
///  --------------------
///  update_x()
///
///     apply correction to solution X:
///        X = X + Y
///     call post_smooth()  [ refresh again first? ]
///
///  --------------------
///  post_smooth()
///
///     apply smoother to A X = B
///     call end_cycle()
///
///  --------------------
///  end_cycle()
///
///     increment iter_
///     if converged exit_solver()
///
///  --------------------
///  exit_solver()
///
///     deallocate temporaries
///     end_compute()
///
///======================================================================
///
/// Required Fields
///
/// - B                          linear system right-hand side
/// - R                          residual R = B - A*X
/// - X                          current solution to A*X = B
/// - Y                          temporary vector for A*X on fine
/// - potential                  computed gravitational potential
/// - density                    density field
/// - acceleration_x             acceleration along X-axis
/// - acceleration_y (rank >= 2) acceleration along Y-axis
/// - acceleration_z (rank >= 3) acceleration along Z-axis

#include "cello.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

//----------------------------------------------------------------------

EnzoMethodGravityMg::EnzoMethodGravityMg 
(FieldDescr * field_descr, int rank,
 double grav_const, int iter_max, double res_tol, int monitor_iter,
 bool is_singular,
 Compute * smooth,
 Restrict * restrict,
 Prolong * prolong,
 int level_min,
 int level_max) 
  : Method(), 
    A_(new EnzoMatrixLaplace),
    smooth_(smooth),
    restrict_(restrict),
    prolong_(prolong),
    is_singular_(is_singular),
    rank_(rank),
    grav_const_(grav_const),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    monitor_iter_(monitor_iter),
    rr_(0),rr0_(0), rr_min_(0),rr_max_(0),
    idensity_(0),  ipotential_(0),
    ib_(0), ir_(0), ix_(0), iy_(0),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    is_active_(0),
    level_(0),
    level_min_(level_min),
    level_max_(level_max)
{
  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  ir_ = field_descr->field_id("R");
  ix_ = field_descr->field_id("X");
  iy_ = field_descr->field_id("Y");
}

//----------------------------------------------------------------------

EnzoMethodGravityMg::~EnzoMethodGravityMg () throw()
{
  delete smooth_;
  delete prolong_;
  delete restrict_;

  smooth_ = NULL;
  prolong_ = NULL;
  restrict_ = NULL;
}


//----------------------------------------------------------------------

void EnzoMethodGravityMg::compute ( Block * block) throw()
{
  set_active(block);

  Field field = block->data()->field();

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  field.size                (&nx_,&ny_,&nz_);
  field.dimensions(idensity_,&mx_,&my_,&mz_);
  field.ghosts    (idensity_,&gx_,&gy_,&gz_);

  int precision = field.precision(idensity_);

  if      (precision == precision_single)    compute_<float>      (enzo_block);
  else if (precision == precision_double)    compute_<double>     (enzo_block);
  else if (precision == precision_quadruple) compute_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityMg()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

extern CkReduction::reducerType r_method_gravity_mg_type;

template <class T>
void EnzoMethodGravityMg::compute_ (EnzoBlock * enzo_block) throw()
//     X = initial guess
//     B = right-hand side
//     R = B - A*X
//     solve(M*Z = R)
//     D = Z
//     shift (B)
{

  set_active(enzo_block);

  iter_ = 0;

  Data * data = enzo_block->data();
  Field field = data->field();

  T * density = (T*) field.values(idensity_);
    
  T * B = (T*) field.values(ib_);
  T * X = (T*) field.values(ix_);
  T * R = (T*) field.values(ir_);
  T * Y = (T*) field.values(iy_);

  // Initialize B, X, R, Y

  if (is_active_) {

    /// X = 0
    /// B = -h^2 * 4 * PI * G * density
    /// R = B ( residual with X = 0 )
    /// Y = 0

    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  B[i] = - 4.0 * (cello::pi) * grav_const_ * density[i];
	  X[i] = 0.0;
	  R[i] = B[i];
	  Y[i] = 0.0;
	}
      }
    }
  }

  pre_smooth(enzo_block,level_max_);
  mg_end<T>(enzo_block,return_unknown);

}


//----------------------------------------------------------------------

void EnzoMethodGravityMg::pre_smooth (EnzoBlock * enzo_block, int level)
{
  smooth_->compute(enzo_block);
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::compute_residual (EnzoBlock * enzo_block, int level)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::restrict_residual (EnzoBlock * enzo_block, int level)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::coarse_solve (EnzoBlock * enzo_block, int level)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::prolong_correction (EnzoBlock * enzo_block, int level)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::update_solution (EnzoBlock * enzo_block, int level)
{
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::post_smooth (EnzoBlock * enzo_block, int level)
{
}

//======================================================================

void EnzoMethodGravityMg::monitor_output_(EnzoBlock * enzo_block) throw()
{

  Monitor * monitor = enzo_block->simulation()->monitor();

  monitor->print ("Enzo", "MG iter %04d  rr %g [%g %g]",
		  iter_,
		  (double)(rr_    / rr0_),
		  (double)(rr_min_/ rr0_),
		  (double)(rr_max_/ rr0_));
}

//----------------------------------------------------------------------

// void EnzoBlock::enzo_matvec_()
// {
//   EnzoMethodGravityMg * method = 
//     static_cast<EnzoMethodGravityMg*> (this->method());

//   Field field = data()->field();
//   int precision = field.precision(field.field_id("density")); assuming 

//   EnzoBlock * enzo_block = static_cast<EnzoBlock*> (this);

//   if      (precision == precision_single)    
//     method->mg_shift_1<float>(enzo_block);
//   else if (precision == precision_double)    
//     method->mg_shift_1<double>(enzo_block);
//   else if (precision == precision_quadruple) 
//     method->mg_shift_1<long double>(enzo_block);
//   else 
//     ERROR1("EnzoMethodGravityMg()", "precision %d not recognized", precision);
// }

//----------------------------------------------------------------------


template <class T>
void EnzoMethodGravityMg::mg_end (EnzoBlock * enzo_block,int retval) throw ()
{
  set_active(enzo_block);

  enzo_block->clear_refresh();

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X         = (T*) field.values(ix_);
  T * potential = (T*) field.values(ipotential_);

  copy_(potential,X,mx_,my_,mz_,is_active_);

  FieldDescr * field_descr = field.field_descr();

  EnzoComputeAcceleration compute_acceleration (field_descr,rank_, true,2);

  compute_acceleration.compute(enzo_block);

  if (enzo_block->index().is_root()) {
    monitor_output_ (enzo_block);
  }

  enzo_block->compute_done();
}

//----------------------------------------------------------------------
