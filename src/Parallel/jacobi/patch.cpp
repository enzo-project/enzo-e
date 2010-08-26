#include <stdio.h>
#include <assert.h>

#include "patch.hpp"
#include "parallel.def"


//----------------------------------------------------------------------

Patch::Patch(int block_count, int block_size, CProxy_Main main_proxy) 
  : main_proxy_ (main_proxy),
    block_count_(block_count),
    block_size_ (block_size),
    count_receive_(6),
    values_(0),
    cycle_values_(0)
{
  allocate_(block_size);
  for (int i=0; i<6; i++) {
    cycle_ghosts_[i] = 0;
  }
  initialize_();
}

//----------------------------------------------------------------------

Patch::~Patch()
{
  delete [] values_;
}

//----------------------------------------------------------------------

Patch::Patch(CkMigrateMessage *) 
{
}

//----------------------------------------------------------------------

void Patch::allocate_(int n) 
{
  block_size_ = n;
  int n1 = n+1;
  values_ = new double [n1*n1*n1];
}

//----------------------------------------------------------------------

void Patch::initialize_()
{
  double xm = 1.0 * thisIndex.x / block_count_;
  double ym = 1.0 * thisIndex.y / block_count_;
  double zm = 1.0 * thisIndex.z / block_count_;
  double xp = 1.0 * (thisIndex.x + 1)/ block_count_;
  double yp = 1.0 * (thisIndex.y + 1) / block_count_;
  double zp = 1.0 * (thisIndex.z + 1) / block_count_;

  for (int iz=1; iz<block_size_-1; iz++) {
    double z = zm + (iz+0.5)/block_size_*(zp-zm);
    for (int iy=1; iy<block_size_-1; iy++) {
      double y = ym + (iy+0.5)/block_size_*(yp-ym);
      for (int ix=1; ix<block_size_-1; ix++) {
	double x = xm + (ix+0.5)/block_size_*(xp-xm);
	int i = ix + block_size_*(iy + block_size_*iz);
	double r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5);
	values_[i] = (r2 < 0.125) ? 1.0 : 0.0;
      }
    }
  }

}

//----------------------------------------------------------------------

void Patch::p_receive(int axis, int face, int n, double * buffer_ghost)
{
  // Update ghosts on given face from buffer

  buffer_to_ghost_(axis, face, buffer_ghost);

  compute_();
}

//----------------------------------------------------------------------
void Patch::compute_()
{
  if (--count_receive_ <= 0) {
    CkPrintf ("Compute %d %d\n",CkMyPe(),cycle_values_);
    ++ cycle_values_;
    count_receive_ = 6;


    main_proxy_.p_next();
  }
  
}

//----------------------------------------------------------------------

void Patch::p_evolve() 
{

  // Update ghost zones

  double * buffer_face[3][2];

  int n3[3]  = { block_size_, block_size_, block_size_ };
  int ng3[3] = { 1,1,1 };

  for (int axis=0; axis<3; axis++) {
    int ng = ng3[axis];
    int n1 = n3[(axis+1)%3];
    int n2 = n3[(axis+2)%3];
    for (int face=0; face<2; face++) {
      buffer_face[axis][face] = new double [ng*n1*n2];
      face_to_buffer_(axis,face,buffer_face[axis][face]);
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

  // (...execution continued in last receive executed locally)
}

//======================================================================

void Patch::face_to_buffer_(int axis, int face, double * buffer)
{
  // TEST BOUNDARY

  // Generalize block size

  int nx = block_size_;
  int ny = block_size_;
  int nz = block_size_;

  // Generalize number of ghost zones

  int ngx = 1;
  int ngy = 1;
  int ngz = 1;

  if (axis == 0 && face == 0) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<ngx; ix++) {
	  int k = ix+ngx;
	  int iv =  k + nx*(iy + ny*iz);
	  int ig = iy + ny*(iz + nz*ix);
	  assert(0 <= ig && ig < ngx*ny*nz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  buffer[ig] = values_[iv];
	}
      }
    }
  } else if (axis == 0 && face == 1) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<ngx; ix++) {
	  int k = ix+nx-2*ngx;
	  int iv =  k + nx*(iy + ny*iz);
	  int ig = iy + ny*(iz + nz*ix);
	  assert(0 <= ig && ig < ngx*ny*nz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  buffer[ig] = values_[iv];
	}
      }
    }
  } else if (axis == 1 && face == 0) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ngy; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int k = iy+ngy;
	  int iv = ix + nx*( k + ny*iz);
	  int ig = iz + nz*(ix + nx*iy);
	  assert(0 <= ig && ig < nx*ngy*nz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  buffer[ig] = values_[iv];
	}
      }
    }
  } else if (axis == 1 && face == 1) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ngy; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int k = iy+ny-2*ngy;
	  int iv = ix + nx*( k + ny*iz);
	  int ig = iz + nz*(ix + nx*iy);
	  assert(0 <= ig && ig < nx*ngy*nz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  buffer[ig] = values_[iv];
	}
      }
    }
  } else if (axis == 2 && face == 0) {
    for (int iz=0; iz<ngz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int k = iz+ngz;
	  int iv = ix + nx*(iy + ny* k);
	  int ig = ix + nx*(iy + ny*iz);
	  assert(0 <= ig && ig < nx*ny*ngz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  buffer[ig] = values_[iv];
	}
      }
    }
  } else if (axis == 2 && face == 1) {
    for (int iz=0; iz<ngz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int k = iz+nz-2*ngz;
	  int iv = ix + nx*(iy + ny* k);
	  int ig = ix + nx*(iy + ny*iz);
	  assert(0 <= ig && ig < nx*ny*ngz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  buffer[ig] = values_[iv];
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void Patch::buffer_to_ghost_(int axis, int face, double * buffer)
{
  // TEST BOUNDARY

  // Generalize block size

  int nx = block_size_;
  int ny = block_size_;
  int nz = block_size_;

  // Generalize number of ghost zones

  int ngx = 1;
  int ngy = 1;
  int ngz = 1;

  if (axis == 0) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<ngx; ix++) {
	  int k = ix + face*(nx-ngx);
	  int iv = k  + nx*(iy + ny*iz);
	  int ig = iy + ny*(iz + nz*ix);
	  assert(0 <= ig && ig < ngx*ny*nz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  values_[iv] = buffer[ig];
	}
      }
    }
  } else if (axis == 1) {
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ngy; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int k = iy + face*(ny-ngy);
	  int iv = ix + nx*( k + ny*iz);
	  int ig = iz + nz*(ix + nx*iy);
	  assert(0 <= ig && ig < nx*ngy*nz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  values_[iv] = buffer[ig];
	}
      }
    }
  } else if (axis == 2) {
    for (int iz=0; iz<ngz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  int k = iz + face*(nz-ngz);
	  int iv = ix + nx*(iy + ny* k);
	  int ig = ix + nx*(iy + ny*iz);
	  assert(0 <= ig && ig < nx*ny*ngz);
	  assert(0 <= iv && iv < nx*ny*nz);
	  values_[iv] = buffer[ig];
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void Patch::print_ ()
{
  double xm = 1.0 * thisIndex.x / block_count_;
  double ym = 1.0 * thisIndex.y / block_count_;
  double zm = 1.0 * thisIndex.z / block_count_;
  double xp = 1.0 * (thisIndex.x + 1)/ block_count_;
  double yp = 1.0 * (thisIndex.y + 1) / block_count_;
  double zp = 1.0 * (thisIndex.z + 1) / block_count_;
  int ix0 = thisIndex.x * block_size_;
  int iy0 = thisIndex.y * block_size_;
  int iz0 = thisIndex.z * block_size_;
  for (int iz=0; iz<block_size_; iz++) {
    double z = zm + (iz+0.5)/block_size_*(zp-zm);
    for (int iy=0; iy<block_size_; iy++) {
      double y = ym + (iy+0.5)/block_size_*(yp-ym);
      for (int ix=0; ix<block_size_; ix++) {
	double x = xm + (ix+0.5)/block_size_*(xp-xm);
	int i = ix + block_size_*(iy + block_size_*iz);
	CkPrintf ("%03d %03d %03d %g\n",ix+ix0,iy+iy0,iz+iz0,values_[i]);
      }
    }
  }
}
