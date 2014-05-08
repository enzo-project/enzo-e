// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHeat.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Declaration of EnzoMethodForwardEuler
///           forward Euler

#ifndef ENZO_ENZO_METHOD_HEAT_HPP
#define ENZO_ENZO_METHOD_HEAT_HPP

class EnzoMethodHeat : public Method {

  /// @class    EnzoMethodHeat
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve heat equation
  /// using forward Euler method

public: // interface

  /// Create a new EnzoMethodHeat object
  EnzoMethodHeat(double alpha, double courant);

  EnzoMethodHeat() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodHeat);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodHeat (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( CommBlock * comm_block) const throw();

  /// Compute maximum timestep for this method
  virtual double timestep ( CommBlock * comm_block) const throw();

  virtual std::string name() const throw() { return "heat"; }

protected: // methods

  template <class T>
  void compute_
  (T * Unew,
   int ndx, int ndy, int ndz,
   int nx,  int ny,  int nz,
   int gx,  int gy,  int gz,
   double dt, double dx, double dy, double dz,
   int dim) const throw();

protected: // attributes

  /// Thermal diffusivity
  double alpha_;

  /// Courant safety number
  double courant_;

};

#endif /* ENZO_ENZO_METHOD_HEAT_HPP */

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

void EnzoMethodHeat::compute
( CommBlock * comm_block) const throw()
{
  Block            * block = comm_block->block();
  FieldBlock * field_block =      block->field_block();
  const FieldDescr * field_descr = comm_block->field_descr();

  const int id_temp   = field_descr->field_id ("temperature");
  void * temp   = field_block->field_values (id_temp);

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  int gx,gy,gz;
  field_descr->ghosts (id_temp,&gx,&gy,&gz);

  int ndx = nx+2*gx;
  int ndy = ny+2*gy;
  int ndz = nz+2*gz;

  double dt = comm_block->dt();

  double dx,dy,dz;
  comm_block->cell_width (&dx,&dy,&dz);

  const int precision = field_descr->precision (id_temp);

  const int dim = comm_block->simulation()->dimension();

  switch (precision) {
  case precision_single:
    compute_ ((float*)  temp, ndx,ndy,ndz, nx,ny,nz, gx,gy,gz, dt, dx,dy,dz, dim);
    break;
  case precision_double:
    compute_ ((double*) temp, ndx,ndy,ndz, nx,ny,nz, gx,gy,gz, dt, dx,dy,dz, dim);
    break;
  default:
    break;
  }
}

//----------------------------------------------------------------------

double EnzoMethodHeat::timestep
(
 CommBlock *        comm_block
 ) const throw()
{
  
  double dx = std::numeric_limits<double>::max();
  double dy = std::numeric_limits<double>::max();
  double dz = std::numeric_limits<double>::max();
  comm_block->cell_width (&dx,&dy,&dz);

  const int dim = comm_block->simulation()->dimension();
  double dm = dx;
  if (dim >= 2) dm = std::min(dm,dy);
  if (dim >= 3) dm = std::min(dm,dz);

  return 0.5*courant_*dm*dm/alpha_;
}

//======================================================================

template <class T>
void EnzoMethodHeat::compute_
(T * Unew, 
 int ndx, int ndy, int ndz,
 int nx,  int ny,  int nz, 
 int gx,  int gy,  int gz, 
 double dt, double dx, double dy, double dz,
 int dim) const throw()
{

  const int idx = 1;
  const int idy = ndx;
  const int idz = ndx*ndy;

  double dxi = 1.0/(dx*dx);
  double dyi = 1.0/(dy*dy);
  double dzi = 1.0/(dz*dz);

  T * U = new T [ndx*ndy*ndz];
  for (int i=0; i<ndx*ndy*ndz; i++) U[i]=Unew[i];

  if (dim == 1) {

    for (int ix=gx; ix<nx+gx; ix++) {

      int i = ix;

      double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);

      Unew[i] = U[i] + alpha_*dt*(Uxx);

    }

  } else if (dim == 2) {

    for (int iy=gy; iy<ny+gy; iy++) {
      for (int ix=gy; ix<nx+gy; ix++) {

	int i = ix + ndx*iy;

	double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	double Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	
	Unew[i] = U[i] + alpha_*dt*(Uxx + Uyy);

      }
    }

  } else if (dim == 3) {

    for (int iz=gz; iz<nz+gz; iz++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int ix=gx; ix<nx+gx; ix++) {

	  int i = ix + ndx*(iy + ndy*iz);

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
