#include <stdio.h>

#include "patch.hpp"
#include "parallel.def"


//----------------------------------------------------------------------

Patch::Patch() 
  : block_count_(0),
    block_size_(0),
    count_receive_(6),
    values_(0),
    cycle_values_(0)
{
  CkPrintf ("Patch()\n");
  for (int i=0; i<6; i++) {
    ghosts_[i]       = 0;
    cycle_ghosts_[i] = 0;
  }
}

//----------------------------------------------------------------------

Patch::~Patch()
{
  CkPrintf ("~Patch()\n");
  delete [] values_;
  for (int i=0; i<6; i++) {
    delete [] ghosts_[i];
  }
}

//----------------------------------------------------------------------

Patch::Patch(CkMigrateMessage *) 
{
  PARALLEL_PRINTF ("Patch(migrate)\n"); 
}

//----------------------------------------------------------------------

void Patch::allocate_(int n) 
{
  block_size_ = n;
  CkPrintf ("allocate_()\n");
  int n1 = n+1;
  values_ = new double [n1*n1*n1];
  for (int i=0; i<6; i++) {
    ghosts_[i] = new double [n1*n1];
  }
}

//----------------------------------------------------------------------

void Patch::initialize_()
{
  CkPrintf ("initialize(%d %d %d)\n",thisIndex.x,thisIndex.y,thisIndex.z);

  double xm = 1.0 * thisIndex.x / block_count_;
  double ym = 1.0 * thisIndex.y / block_count_;
  double zm = 1.0 * thisIndex.z / block_count_;
  double xp = 1.0 * (thisIndex.x + 1)/ block_count_;
  double yp = 1.0 * (thisIndex.y + 1) / block_count_;
  double zp = 1.0 * (thisIndex.z + 1) / block_count_;
  CkPrintf ("initialize(%g %g  %g %g  %g %g)\n",
	    xm,xp,ym,yp,zm,zp);

  for (int iz=1; iz<block_size_-1; iz++) {
    double z = zm + (iz+0.5)/block_size_*(zp-zm);
    for (int iy=1; iy<block_size_-1; iy++) {
      double y = ym + (iy+0.5)/block_size_*(yp-ym);
      for (int ix=1; ix<block_size_-1; ix++) {
	double x = xm + (ix+0.5)/block_size_*(xp-xm);
	int i = ix + block_size_*(iy + block_size_*iz);
	values_[i] = initial_(x,y,z);
      }
    }
  }

}

//----------------------------------------------------------------------

double Patch::initial_(double x, double y, double z)
{
  double r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5);
  return (r2 < 0.125) ? 1.0 : 0.0;
}

//----------------------------------------------------------------------

void Patch::p_receive(int face, int dir, int n, double * buffer_ghost)
{
  CkPrintf ("p_receive(%d)\n",count_receive_);

  // Update ghosts on given face from buffer

  buffer_to_ghost_(face, dir, buffer_ghost);
  
  // check for all ghosts received
  if (--count_receive_ <= 0) {
    print_();
  }
}

//----------------------------------------------------------------------

void Patch::p_evolve(int block_count, int block_size) 
{

  CkPrintf ("p_evolve()\n");

  block_size_ = block_size;
  block_count_ = block_count;

  allocate_(block_size);

  initialize_();

  double * buffer_face[3][2];

  for (int axis=0; axis<3; axis++) {
    for (int dir=0; dir<2; dir++) {
      buffer_face[axis][dir] = new double [block_size*block_size];
      face_to_buffer_(axis,dir,buffer_face[axis][dir]);
    }
  }

  CProxy_Patch patch = thisProxy;

  int ix = thisIndex.x;
  int iy = thisIndex.y;
  int iz = thisIndex.z;

  int ixm = (ix - 1 + block_count_) % block_count_;
  int iym = (iy - 1 + block_count_) % block_count_;
  int izm = (iz - 1 + block_count_) % block_count_;
  int ixp = (ix + 1) % block_count_;
  int iyp = (iy + 1) % block_count_;
  int izp = (iz + 1) % block_count_;

  patch(ixp,iy,iz).p_receive(0,0,block_size_,buffer_face[0][1]);
  patch(ixm,iy,iz).p_receive(0,1,block_size_,buffer_face[0][0]);
  patch(ix,iyp,iz).p_receive(1,0,block_size_,buffer_face[1][1]);
  patch(ix,iym,iz).p_receive(1,1,block_size_,buffer_face[1][0]);
  patch(ix,iy,izp).p_receive(2,0,block_size_,buffer_face[2][1]);
  patch(ix,iy,izm).p_receive(2,1,block_size_,buffer_face[2][0]);

  // ...execution continued in last receive executed locally
}

//======================================================================

void Patch::face_to_buffer_(int axis, int dir, double * buffer_face)
{
  CkPrintf ("face_to_buffer_()\n");

#error  Implementing single-loop face_to_buffer (3D to 2D)

  if (axis == 0 && dir == 0) ix = 1;
  if (axis == 1 && dir == 0) iy = 1;
  if (axis == 2 && dir == 0) iz = 1;
  if (axis == 0 && dir == 1) ix = block_size_-2;
  if (axis == 1 && dir == 1) iy = block_size_-2;
  if (axis == 2 && dir == 1) iz = block_size_-2;

  int d1,d2;
  d1 = (axis == 0) ? 0 : 1;
  d2 = (axis == 0) ? 0 : block_size_;

  for (int i2=0; i2<block_size_; i2++) {
    for (int i1=0; i1<block_size_; i1++) {
      int i_buffer = i1 + block_size_*i2;
      int i_values = i1*dx + i2*dy;
      buffer_face[i_buffer] = values_[i_values];
    }
  }
}

//----------------------------------------------------------------------

void Patch::buffer_to_ghost_(int face, double * buffer_ghost)
{
  CkPrintf ("buffer_to_ghost_()\n");
  int ix,iy,iz;

  // Exit if boundary

  // Initialize constant index

  switch (face) {
  case face_xm: ix = 0; break;
  case face_ym: iy = 0; break;
  case face_zm: iz = 0; break;
  case face_xp: ix = block_size_-1; break;
  case face_yp: iy = block_size_-1; break;
  case face_zp: iz = block_size_-1; break;
  }

  // Copy buffer to ghosts

  switch (face) {
  case face_xm:
  case face_xp:
    for (int iz=0; iz<block_size_; iz++) {
      for (int iy=0; iy<block_size_; iy++) {
	int i_buffer = iy + block_size_*iz;
	int i_values = ix + block_size_*(iy + block_size_*iz);
	values_[i_values] = buffer_ghost[i_buffer];
      }
    }
    break;
  case face_ym:
  case face_yp:
    for (int iz=0; iz<block_size_; iz++) {
      for (int ix=0; ix<block_size_; ix++) {
	int i_buffer = iz + block_size_*ix;
	int i_values = ix + block_size_*(iy + block_size_*iz);
	values_[i_values] = buffer_ghost[i_buffer];
      }
    }
    break;
  case face_zm:
  case face_zp:
    for (int iy=0; iy<block_size_; iy++) {
      for (int ix=0; ix<block_size_; ix++) {
	int i_buffer = ix + block_size_*iy;
	int i_values = ix + block_size_*(iy + block_size_*iz);
	values_[i_values] = buffer_ghost[i_buffer];
      }
    }
    break;
  }
}

//----------------------------------------------------------------------

void Patch::print_ ()
{
  CkPrintf ("print_()\n");
  double xm = 1.0 * thisIndex.x / block_count_;
  double ym = 1.0 * thisIndex.y / block_count_;
  double zm = 1.0 * thisIndex.z / block_count_;
  double xp = 1.0 * (thisIndex.x + 1)/ block_count_;
  double yp = 1.0 * (thisIndex.y + 1) / block_count_;
  double zp = 1.0 * (thisIndex.z + 1) / block_count_;
  for (int iz=0; iz<block_size_; iz++) {
    double z = zm + (iz+0.5)/block_size_*(zp-zm);
    for (int iy=0; iy<block_size_; iy++) {
      double y = ym + (iy+0.5)/block_size_*(yp-ym);
      for (int ix=0; ix<block_size_; ix++) {
	double x = xm + (ix+0.5)/block_size_*(xp-xm);
	int i = ix + block_size_*(iy + block_size_*iz);
	CkPrintf ("%g %g %g %g\n",x,y,z,values_[i]);
      }
    }
  }
}
