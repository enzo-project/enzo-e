#include <stdio.h>

#include "patch.hpp"
#include "parallel.def"


//----------------------------------------------------------------------

Patch::Patch() 
  : n_(0),
    values_(0),
    cycle_values_(0),
    cycle_ghosts_(0)
{
  CkPrintf ("PATCH()\n");
  for (int i=0; i<6; i++) {
    ghosts_[i] = 0;
  }
}

//----------------------------------------------------------------------

Patch::~Patch()
{
  delete [] values_;
  for (int i=0; i<6; i++) {
    delete [] ghosts_[i];
  }
}

//----------------------------------------------------------------------

Patch::Patch(CkMigrateMessage *) 
{
  PARALLEL_PRINTF ("PATCH(MIGRATE)\n"); 
}

//----------------------------------------------------------------------

void Patch::allocate(int n) 
{
  n_ = n;
  CkPrintf ("ALLOCATE()\n");
  int n1 = n+1;
  values_ = new double [n1*n1*n1];
  for (int i=0; i<6; i++) {
    ghosts_[i] = new double [n1*n1];
  }
}

//----------------------------------------------------------------------

void Patch::neighbors (CProxy_Patch xm, CProxy_Patch xp,
		       CProxy_Patch ym, CProxy_Patch yp,
		       CProxy_Patch zm, CProxy_Patch zp)
{
  neighbor_[neighbor_xm] = xm;
  neighbor_[neighbor_xp] = xp;
  neighbor_[neighbor_ym] = ym;
  neighbor_[neighbor_yp] = yp;
  neighbor_[neighbor_zm] = zm;
  neighbor_[neighbor_zp] = zp;
}

//----------------------------------------------------------------------

void Patch::initialize()
{
  CkPrintf ("INITIALIZE()\n");
  printf ("%d %d %d\n",thisIndex.x,thisIndex.y,thisIndex.z);
//   for (int iz=0; iz<n_; iz++) {
//     double z = zm + (iz+0.5)/n_*(zp-zm);
//     for (int iy=0; iy<n_; iy++) {
//       double y = ym + (iy+0.5)/n_*(yp-ym);
//       for (int ix=0; ix<n_; ix++) {
// 	double x = xm + (ix+0.5)/n_*(xp-xm);
// 	int i = ix + n_*(iy + n_*iz);
// 	values_[i] = initial_(x,y,z);
//       }
//     }
//   }
}

//----------------------------------------------------------------------

double Patch::initial_(double x, double y, double z)
{
}

//----------------------------------------------------------------------

void Patch::receive(int type, int n, double * buffer)
{
  switch (type) {
  case neighbor_xm:  CkPrintf ("RECEIVE XM\n"); break;
  case neighbor_xp:  CkPrintf ("RECEIVE XP\n"); break;
  case neighbor_ym:  CkPrintf ("RECEIVE YM\n"); break;
  case neighbor_yp:  CkPrintf ("RECEIVE YP\n"); break;
  case neighbor_zm:  CkPrintf ("RECEIVE ZM\n"); break;
  case neighbor_zp:  CkPrintf ("RECEIVE ZP\n"); break;
  }
}

//----------------------------------------------------------------------

void Patch::advance(int n) 
{

  CkPrintf ("ADVANCE()\n");

  allocate(n);

  initialize();

  double * buffer = new double [n*n];

  prepare_(neighbor_xp,buffer);

//   for (int iz=0; iz<nb; iz++) {
//     for (int iy=0; iy<nb; iy++) {
//       for (int ix=0; ix<nb; ix++) {
// 	if (ix<nb-1) patch(ix+1,iy,iz).receive(neighbor_xp,n,buffer);
// 	if (ix>0)    patch(ix-1,iy,iz).receive(neighbor_xm,n,buffer);
// 	if (iy<nb-1) patch(ix,iy+1,iz).receive(neighbor_yp,n,buffer);
// 	if (iy>0)    patch(ix,iy-1,iz).receive(neighbor_ym,n,buffer);
// 	if (iz<nb-1) patch(ix,iy,iz+1).receive(neighbor_zp,n,buffer);
// 	if (iz>0)    patch(ix,iy,iz-1).receive(neighbor_zm,n,buffer);
//       }
//     }
//   }

  unpack_(neighbor_xp,buffer);



}

//======================================================================

void Patch::prepare_(int type, double * buffer)
{
  CkPrintf ("PREPARE\n");
}

//----------------------------------------------------------------------

void Patch::unpack_(int type, double * buffer)
{
  CkPrintf ("UNPACK\n");
}

