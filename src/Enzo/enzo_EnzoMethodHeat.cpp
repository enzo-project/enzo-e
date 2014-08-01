// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHeat.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodHeat class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodHeat::EnzoMethodHeat (double alpha, double courant) 
  : Method(), alpha_(alpha), courant_(courant)
{
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

void EnzoMethodHeat::compute ( CommBlock * comm_block) throw()
{

  if (!comm_block->is_leaf()) return;

  initialize_(comm_block);

  const int id = field_id ("temperature");
  void *     t = field_array (id);
  const int  p = field_precision (id);

  if      (p == precision_single)    compute_ ((float*)  t);
  else if (p == precision_double)    compute_ ((double*) t);
  else if (p == precision_quadruple) compute_ ((long double*) t);
  else 
    ERROR1("EnzoMethodHeat()", "precision %d not recognized", p);
}

//----------------------------------------------------------------------

double EnzoMethodHeat::timestep ( CommBlock * comm_block ) throw()
{
  initialize_(comm_block);

  const int rank = this->rank();

  double hx,hy,hz;
  cell_width (&hx,&hy,&hz);

  double h_min = std::numeric_limits<double>::max();
  if (rank >= 1) h_min = std::min(h_min,hx);
  if (rank >= 2) h_min = std::min(h_min,hy);
  if (rank >= 3) h_min = std::min(h_min,hz);

  return 0.5*courant_*h_min*h_min/alpha_;
}

//======================================================================

template <class T>
void EnzoMethodHeat::compute_ (T * Unew) const throw()
{

  const int id = field_id ("temperature");

  int gx,gy,gz;
  ghost_depth (id,&gx,&gy,&gz);

  int mx,my,mz;
  array_dimension (id,&mx,&my,&mz);

  int nx,ny,nz;
  block_size (&nx,&ny,&nz);

  // Initialize array increments
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my;

  // Precompute ratios
  double hx,hy,hz;
  cell_width (&hx,&hy,&hz);
  double dxi = 1.0/(hx*hx);
  double dyi = 1.0/(hy*hy);
  double dzi = 1.0/(hz*hz);

  const int m = mx*my*mz;

  const int rank = this->rank();

  const double dt = time_step();

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
