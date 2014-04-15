// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHeat.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodHeat class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodHeat::EnzoMethodHeat (double alpha) 
  : Method(),
    alpha_(alpha)
{
  printf ("%s:%d Created EnzoMethodHeat\n",__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoMethodHeat::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodHeat::compute_block
(
 FieldDescr * field_descr, CommBlock * comm_block) throw()
{
  Block            * block = comm_block->block();
  FieldBlock * field_block =      block->field_block();

  const int id_temp   = field_descr->field_id ("temperature");

  void * temp   = field_block->field_unknowns (id_temp);

  int ndx,ndy,ndz;
  field_block->size(&ndx,&ndy,&ndz);

  int ngx,ngy,ngz;
  field_descr->ghosts (id_temp,&ngx,&ngy,&ngz);

  int nx = ndx-2*ngx;
  int ny = ndy-2*ngy;
  int nz = ndz-2*ngz;

  double dt,dx,dy,dz;
  comm_block->cell_width (&dx,&dy,&dz);

  const int precision = field_descr->precision (id_temp);

  switch (precision) {
  case precision_single:
    compute ((float*)  temp, ndx,ndy,ndz, nx,ny,nz, dt, dx,dy,dz);
    break;
  case precision_double:
    compute ((double*) temp, ndx,ndy,ndz, nx,ny,nz, dt, dx,dy,dz);
    break;
  default:
    break;
  }
}

template <class T>
void EnzoMethodHeat::compute 
(T * Unew, 
 int ndx, int ndy, int ndz,
 int nx,  int ny,  int nz, 
 double dt, double dx, double dy, double dz) const throw()
{
  {
    const int idx = 1;
    const int idy = ndx;
    const int idz = ndx*ndy;
    double dxi = 1.0/dx;
    double dyi = 1.0/dy;
    double dzi = 1.0/dz;
    T * U = new T [ndx*ndy*ndz];
    for (int i=0; i<ndx*ndy*ndz; i++) U[i]=Unew[i];
    
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i = ix + ndx*(iy + ndy*iz);
	  double Uxx = dxi*(U[i-idx] - 2*U[i] + U[i+idx]);
	  double Uyy = dyi*(U[i-idy] - 2*U[i] + U[i+idy]);
	  double Uzz = dzi*(U[i-idz] - 2*U[i] + U[i+idz]);
	  Unew[i] = U[i] + alpha_*dt*(Uxx + Uyy + Uzz);
	}
      }
    }
    delete [] U;
  }
}
