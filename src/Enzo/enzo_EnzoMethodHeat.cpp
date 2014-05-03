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

  initialize_(comm_block);

  const int id = field_id_["temperature"];
  void *     t = field_array_[id];
  const int  p = field_precision_[id];

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
  double h = hx_;
  if (rank_ >= 2) h = std::min(h,hy_);
  if (rank_ >= 3) h = std::min(h,hz_);

  return 0.5*courant_*h*h/alpha_;
}

//======================================================================

template <class T>
void EnzoMethodHeat::compute_ (T * Unew) const throw()
{

  const int id = field_id_.at("temperature");

  const int gx = gx_[id];
  const int gy = gy_[id];
  const int gz = gz_[id];

  const int mx = mx_[id];
  const int my = my_[id];
  const int mz = mz_[id];

  // Initialize array increments
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my;

  // Precompute ratios
  double dxi = 1.0/(hx_*hx_);
  double dyi = 1.0/(hy_*hy_);
  double dzi = 1.0/(hz_*hz_);

  const int m = mx*my*mz;
  T * U = new T [m];
  for (int i=0; i<m; i++) U[i]=Unew[i];

  if (rank_ == 1) {

    for (int ix=gx; ix<nx_+gx; ix++) {

      int i = ix;

      double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);

      Unew[i] = U[i] + alpha_*dt_*(Uxx);

    }

  } else if (rank_ == 2) {

    for (int iy=gy; iy<ny_+gy; iy++) {
      for (int ix=gy; ix<nx_+gy; ix++) {

	int i = ix + mx*iy;

	double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	double Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	
	Unew[i] = U[i] + alpha_*dt_*(Uxx + Uyy);

      }
    }

  } else if (rank_ == 3) {

    for (int iz=gz; iz<nz_+gz; iz++) {
      for (int iy=gy; iy<ny_+gy; iy++) {
	for (int ix=gx; ix<nx_+gx; ix++) {

	  int i = ix + mx*(iy + my*iz);

	  double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	  double Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	  double Uzz = dzi*(U[i-idz] - 2*U[i] + U[i+idz]);

	  Unew[i] = U[i] + alpha_*dt_*(Uxx + Uyy + Uzz);

	}
      }
    }
  }
  delete [] U;

}
