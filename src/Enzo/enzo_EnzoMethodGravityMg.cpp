// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityMg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityMg class

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
 bool is_singular) 
  : Method(), 
    A_(new EnzoMatrixLaplace),
    is_singular_(is_singular),
    rank_(rank),
    grav_const_(grav_const),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    monitor_iter_(monitor_iter),
    rr0_(0),rr_(0), rr_min_(0),rr_max_(0),
    idensity_(0),  ipotential_(0),
    ib_(0), ix_(0), ir_(0), ic_(0),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    level_(0),
    is_active_(0)
{
  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  ix_ = field_descr->field_id("X");
  ir_ = field_descr->field_id("R");
  ic_ = field_descr->field_id("C");
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
  T * C = (T*) field.values(ic_);

  // Initialize B, X, R, C

  if (is_active_) {

    /// X = 0
    /// B = -h^2 * 4 * PI * G * density
    /// R = B ( residual with X = 0 )
    /// C = 0

    for (int iz=0; iz<mz_; iz++) {
      for (int iy=0; iy<my_; iy++) {
	for (int ix=0; ix<mx_; ix++) {
	  int i = ix + mx_*(iy + my_*iz);
	  B[i] = - 4.0 * (cello::pi) * grav_const_ * density[i];
	  X[i] = 0.0;
	  R[i] = B[i];
	  C[i] = 0.0;
	}
      }
    }
  }

  mg_end<T>(enzo_block,return_unknown);

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
///    if (return == return_converged) {
///       potential = X
///       ==> mg_exit()
///    } else {
///       ERROR (return-)
///    }
{
  set_active(enzo_block);

  enzo_block->clear_refresh();

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X         = (T*) field.values(ix_);
  T * potential = (T*) field.values(ipotential_);

  copy_(potential,X,mx_,my_,mz_,is_active_);

  bool symmetric;
  int order;
  EnzoComputeAcceleration compute_acceleration (field.field_descr(),
						rank_, symmetric = true,
						order=2);

  compute_acceleration.compute(enzo_block);

  if (enzo_block->index().is_root()) {
    monitor_output_ (enzo_block);
  }

  enzo_block->compute_done();
}

//----------------------------------------------------------------------
