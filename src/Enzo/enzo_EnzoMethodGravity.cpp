// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class


#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_METHOD

#ifdef DEBUG_METHOD
#   define TRACE_METHOD(METHOD,BLOCK)					\
  CkPrintf ("%d %s:%d %s TRACE %s %p\n",CkMyPe(),__FILE__,__LINE__, \
	    BLOCK->name().c_str(),METHOD,this);			    \
  fflush(stdout);
#else
#   define TRACE_METHOD(METHOD,BLOCK) /*  */ 
#endif
//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity
(const FieldDescr * field_descr,
 int index_solver,
 double grav_const)
  : Method(),
    index_solver_(index_solver),
    grav_const_(grav_const)
{
  const int num_fields = field_descr->field_count();
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields(num_fields);
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{

  TRACE_METHOD("compute()",block);

  // Initialize the linear system

  Field field = block->data()->field();

  Matrix * A = new EnzoMatrixLaplace;
  const int ix = field.field_id ("potential");

  /// access problem-defining fields for eventual RHS and solution
  const int id  = field.field_id("density");
  const int idt = field.field_id("density_total");
  const int idensity = (idt != -1) ? idt : id;
  const int ib = field.field_id ("B");
  
  // Solve the linear system
  field.scale(ib, -4.0 * (cello::pi) * grav_const_, idensity);

  Solver * solver = block->simulation()->problem()->solver(index_solver_);
  
  // May exit before solve is done...
  solver->set_callback (CkIndex_EnzoBlock::r_method_gravity_continue());

  solver->apply (A, ix, ib, block);

}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_continue()
{

  TRACE_METHOD("r_method_gravity_end()",this);

  // So do refresh with barrier synch (note barrier instead of
  // neighbor synchronization otherwise will conflict with Method
  // refresh ("Charm++ fatal error: mis-matched client callbacks in
  // reduction messages")

  Refresh refresh (4,0,neighbor_leaf, sync_barrier);
  refresh.set_active(is_leaf());
  refresh.add_all_fields(data()->field().field_count());

  refresh_enter(CkIndex_EnzoBlock::r_method_gravity_end(NULL),&refresh);

}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_end(CkReductionMsg * msg)
{
  TRACE_METHOD("r_method_gravity_end()",this);
  
  delete msg;
  
  // BUG: acceleration computed before Solver completes
  
  /// compute acceleration fields from potential
  int order;
  EnzoComputeAcceleration compute_acceleration(data()->field().field_descr(),
					       rank(), order=4);

  compute_acceleration.compute(this);

  // wait for all Blocks before continuing
  compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodGravity::timestep (Block * block) const throw()
{
  Field field = block->data()->field();
  
  int precision = field.precision(0);

  if      (precision == precision_single)
    return timestep_<float>      (block);
  else if (precision == precision_double)
    return timestep_<double>     (block);
  else if (precision == precision_quadruple)
    return timestep_<long double>(block);
  else 
    ERROR1("EnzoMethodGravity()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
double EnzoMethodGravity::timestep_ (Block * block) const throw()
{
  Field field = block->data()->field();

  int nx,ny,nz;
  int mx,my,mz;
  int gx,gy,gz;
  field.size         (&nx,&ny,&nz);
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  T * ax = (T*) field.values ("acceleration_x");
  T * ay = (T*) field.values ("acceleration_y");
  T * az = (T*) field.values ("acceleration_z");

  double dt = std::numeric_limits<double>::max();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);
  
  if (ax) {
    for (int ix=gx; ix<nx+gx; ix++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int iz=gz; iz<nz+gz; iz++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(dt,sqrt(hx/(fabs(ax[i]+1e-20))));
	}
      }
    }
  }
  if (ay) {
    for (int ix=gx; ix<nx+gx; ix++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int iz=gz; iz<nz+gz; iz++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(dt,sqrt(hy/(fabs(ay[i]+1e-20))));
	}
      }
    }
  }
  if (az) {
    for (int ix=gx; ix<nx+gx; ix++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int iz=gz; iz<nz+gz; iz++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(dt,sqrt(hz/(fabs(az[i]+1e-20))));
	}
      }
    }
  }

  return 0.5*dt;
}
