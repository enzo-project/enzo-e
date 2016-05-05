// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmUpdate.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPmUpdate class

#include "cello.hpp"
#include "enzo.hpp"

#define DEBUG_PM_UPDATE

#ifdef DEBUG_PM_UPDATE
#  define TRACE_PM(MESSAGE)						\
  CkPrintf ("%s:%d %s\n",						\
	    __FILE__,__LINE__,MESSAGE);				
#else
#  define TRACE_PM(MESSAGE) /* ... */
#endif

//----------------------------------------------------------------------

EnzoMethodPmUpdate::EnzoMethodPmUpdate ( const FieldDescr * field_descr ) 
  : Method()
{
  TRACE_PM("EnzoMethodPmUpdate()");
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields(field_descr->field_count());

  // PM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPmUpdate::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
}

//----------------------------------------------------------------------

void EnzoMethodPmUpdate::compute ( Block * block) throw()
{
  TRACE_PM("compute()");
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  if (block->is_leaf()) {

    Particle particle = block->data()->particle();
    Field field = block->data()->field();

    int mx,my,mz;
    field.dimensions(0,&mx,&my,&mz);
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);

    const int it = particle.type_index("dark");

    const int rank = block->rank();

    if (rank == 1) {
    } else if (rank == 2) {
      
      double * x = (double *)field.values("velocity_x");
      double * y = (double *)field.values("velocity_y");
      double * vx = (double *)field.values("velocity_x");
      double * vy = (double *)field.values("velocity_y");
      double * ax = (double *)field.values("acceleration_x");
      double * ay = (double *)field.values("acceleration_y");


      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int i = ix + mx*iy;
	  //	  de_t[i]  = de[i];
	}
      }
    } else if (rank == 3) {
    }
  /* 1) v(n) --> v(n+1/2) with a(n+1/2) */
  /* 2) x(n) --> x(n+1) with v(n+1/2)
     /*  if (ProblemType != 23) */
    /*     Grid->UpdateParticlePosition(dt) */
  /* 3) v(n+1/2) --> v(n+1) with a(n+1/2) */

  }

  block->compute_done(); 
  
}

//----------------------------------------------------------------------

double EnzoMethodPmUpdate::timestep ( Block * block ) const throw()
{
  TRACE_PM("timestep()");

  double dt = std::numeric_limits<double>::max();

  return dt;
}
