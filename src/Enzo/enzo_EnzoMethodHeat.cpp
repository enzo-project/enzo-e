// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHeat.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodHeat class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodHeat::EnzoMethodHeat (const FieldDescr * field_descr, double alpha, double courant) 
  : Method(),
    alpha_(alpha),
    courant_(courant)
{
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);

  refresh(ir)->add_all_fields(field_descr->field_count());
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

  const int id_temp = field.field_id ("temperature");

  void *     t = field.values (id_temp);
  const int  p = field.precision (id_temp);

  if      (p == precision_single)    compute_ (block,(float *)t);
  else if (p == precision_double)    compute_ (block,(double*)t);
  else if (p == precision_quadruple) compute_ (block,(long double*) t);
  else 
    ERROR1("EnzoMethodHeat()", "precision %d not recognized", p);
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

  double xm,ym,zm;
  data->lower(&xm,&ym,&zm);
  double xp,yp,zp;
  data->upper(&xp,&yp,&zp);
  double hx,hy,hz;
  field.cell_width(xm,xp,&hx,
		   ym,yp,&hy,
		   zm,zp,&hz);

  double h_min = std::numeric_limits<double>::max();
  if (rank >= 1) h_min = std::min(h_min,hx);
  if (rank >= 2) h_min = std::min(h_min,hy);
  if (rank >= 3) h_min = std::min(h_min,hz);

  return 0.5*courant_*h_min*h_min/alpha_;
}

//======================================================================

template <class T>
void EnzoMethodHeat::compute_ (Block * block,T * Unew) const throw()
{
  Data * data = block->data();
  Field field   =      data->field();

  const int id_temp_ = field.field_id ("temperature");

  int gx,gy,gz;
  field.ghost_depth (id_temp_,&gx,&gy,&gz);

  int mx,my,mz;
  field.dimensions (id_temp_,&mx,&my,&mz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  // Initialize array increments
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my;

  // Precompute ratios dxi,dyi,dzi

  double xm,ym,zm;
  data->lower(&xm,&ym,&zm);
  double xp,yp,zp;
  data->upper(&xp,&yp,&zp);
  double hx,hy,hz;
  field.cell_width(xm,xp,&hx,
		   ym,yp,&hy,
		   zm,zp,&hz);

  double dxi = 1.0/(hx*hx);
  double dyi = 1.0/(hy*hy);
  double dzi = 1.0/(hz*hz);

  const int m = mx*my*mz;

  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  const double dt = timestep(block);

  T * U = new T [m];
  for (int i=0; i<m; i++) U[i]=Unew[i];

  if (rank == 1) {

    for (int ix=gx; ix<nx+gx; ix++) {

      int i = ix;

      double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);

      Unew[i] = U[i] + alpha_*dt*(Uxx);

    }

  } else if (rank == 2) {

    for (int iy=gy; iy<ny+gy; iy++) {
      for (int ix=gy; ix<nx+gy; ix++) {

	int i = ix + mx*iy;

	double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	double Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	
	Unew[i] = U[i] + alpha_*dt*(Uxx + Uyy);

      }
    }

  } else if (rank == 3) {

    for (int iz=gz; iz<nz+gz; iz++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int ix=gx; ix<nx+gx; ix++) {

	  int i = ix + mx*(iy + my*iz);

	  double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	  double Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	  double Uzz = dzi*(U[i-idz] - 2*U[i] + U[i+idz]);

	  Unew[i] = U[i] + alpha_*dt*(Uxx + Uyy + Uzz);

	}
      }
    }
  }

  delete [] U;

}
