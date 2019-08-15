// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHeat.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodHeat class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodHeat::EnzoMethodHeat (double alpha, double courant) 
  : Method(),
    alpha_(alpha),
    courant_(courant)
{
  // Initialize default Refresh object

#ifdef NEW_REFRESH
#else  
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
			     enzo_sync_id_method_heat);
  refresh(ir)->add_field("temperature");
#endif  

}

//----------------------------------------------------------------------

void EnzoMethodHeat::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | alpha_;
  p | courant_;
}

//----------------------------------------------------------------------

void EnzoMethodHeat::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    Field field = block->data()->field();

    enzo_float * T = (enzo_float *) field.values ("temperature");

    compute_ (block,T);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodHeat::timestep ( Block * block ) const throw()
{
  // initialize_(block);

  Data * data = block->data();
  Field field = data->field();

  const int id_temp = field.field_id ("temperature");

  int mx,my,mz;
  field.dimensions (id_temp,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  double h_min = std::numeric_limits<double>::max();
  if (rank >= 1) h_min = std::min(h_min,hx);
  if (rank >= 2) h_min = std::min(h_min,hy);
  if (rank >= 3) h_min = std::min(h_min,hz);

  return 0.5*courant_*h_min*h_min/alpha_;
}

//======================================================================

void EnzoMethodHeat::compute_ (Block * block,enzo_float * Unew) const throw()
{
  Data * data = block->data();
  Field field   =      data->field();

  const int id_temp_ = field.field_id ("temperature");

  int mx,my,mz;
  int gx,gy,gz;

  field.dimensions  (id_temp_,&mx,&my,&mz);
  field.ghost_depth (id_temp_,&gx,&gy,&gz);

  // Initialize array increments
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my;

  // Precompute ratios dxi,dyi,dzi

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  double dxi = 1.0/(hx*hx);
  double dyi = 1.0/(hy*hy);
  double dzi = 1.0/(hz*hz);

  const int m = mx*my*mz;

  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  const double dt = timestep(block);

  enzo_float * U = new enzo_float [m];
  for (int i=0; i<m; i++) U[i]=Unew[i];

  if (rank == 1) {

    for (int ix=gx; ix<mx-gx; ix++) {

      int i = ix;

      enzo_float Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);

      Unew[i] = U[i] + alpha_*dt*(Uxx);

    }

  } else if (rank == 2) {

    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gy; ix<mx-gy; ix++) {

	int i = ix + mx*iy;

	enzo_float Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	enzo_float Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	
	Unew[i] = U[i] + alpha_*dt*(Uxx + Uyy);

      }
    }

  } else if (rank == 3) {

    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {

	  int i = ix + mx*(iy + my*iz);

	  enzo_float Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	  enzo_float Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	  enzo_float Uzz = dzi*(U[i-idz] - 2*U[i] + U[i+idz]);

	  Unew[i] = U[i] + alpha_*dt*(Uxx + Uyy + Uzz);

	}
      }
    }
  }

  delete [] U;

}
