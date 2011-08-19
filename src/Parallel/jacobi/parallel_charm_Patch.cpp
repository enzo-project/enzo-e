
#include <stdio.h>
#include <assert.h>
#include <hdf5.h>

#include "parallel_jacobi.hpp"

#include "parallel.def"

//----------------------------------------------------------------------

CharmPatch::CharmPatch(int num_blocks, int block_size, int cycle_max,
		       CProxy_Main main_proxy) 
  : main_proxy_ (main_proxy),
    values_(0),
    cycle_write_(1),
    cycle_values_(0),
    cycle_max_(cycle_max),
    receives_(6)
{

  // Number of blocks along each axis

  nbx_ = num_blocks;
  nby_ = num_blocks;
  nbz_ = num_blocks;

  // Depth of ghost zone layer along each block face

  ngx_ = 1;
  ngy_ = 1;
  ngz_ = 1;

  // Number of values along each block axis

  nvx_ = block_size;
  nvy_ = block_size;
  nvz_ = block_size;

  // Allocated array size along each block axis

  nax_ = nvx_ + 2*ngx_;
  nay_ = nvy_ + 2*ngy_;
  naz_ = nvz_ + 2*ngz_;
  
  // Lower indices for block values

  ixl_ = ngx_;
  iyl_ = ngy_;
  izl_ = ngz_;

  // Upper indices for block values

  ixu_ = ngx_ + nvx_;
  iyu_ = ngy_ + nvy_;
  izu_ = ngz_ + nvz_;

  allocate_();

  initialize_();
}

//----------------------------------------------------------------------

CharmPatch::~CharmPatch()
{
  delete [] values_;
  delete [] buffer_[0][0];
  delete [] buffer_[0][1];
  delete [] buffer_[1][0];
  delete [] buffer_[1][1];
  delete [] buffer_[2][0];
  delete [] buffer_[2][1];
}

//----------------------------------------------------------------------

CharmPatch::CharmPatch(CkMigrateMessage *) 
  : receives_(6)
{
}

//----------------------------------------------------------------------

void CharmPatch::allocate_() 
{
  values_       = new double [nax_*nay_*naz_];
  buffer_[0][0] = new double [ngx_*nay_*naz_];
  buffer_[0][1] = new double [ngx_*nay_*naz_];
  buffer_[1][0] = new double [nax_*ngy_*naz_];
  buffer_[1][1] = new double [nax_*ngy_*naz_];
  buffer_[2][0] = new double [nax_*nay_*ngz_];
  buffer_[2][1] = new double [nax_*nay_*ngz_];
}

// //----------------------------------------------------------------------

// int CharmPatch::id_()
// {
//   return thisIndex.x + nbx_*(thisIndex.y + nby_*thisIndex.z); 
// };

//----------------------------------------------------------------------

void CharmPatch::initialize_()
{
  const double R2 = 0.25*0.25;
  double xm = 1.0 * thisIndex.x / nbx_;
  double ym = 1.0 * thisIndex.y / nby_;
  double zm = 1.0 * thisIndex.z / nbz_;
  double xp = 1.0 * (thisIndex.x + 1) / nbx_;
  double yp = 1.0 * (thisIndex.y + 1) / nby_;
  double zp = 1.0 * (thisIndex.z + 1) / nbz_;

  for (int iz=izl_; iz<izu_; iz++) {

    double z  = zm + (iz-izl_+0.5)/nvz_*(zp-zm);
    double z2 = (z-0.5)*(z-0.5);

    for (int iy=iyl_; iy<iyu_; iy++) {

      double y  = ym + (iy-iyl_+0.5)/nvy_*(yp-ym);
      double y2 = (y-0.5)*(y-0.5);

      for (int ix=ixl_; ix<ixu_; ix++) {

	double x  = xm + (ix-ixl_+0.5)/nvx_*(xp-xm);
	double x2 = (x-0.5)*(x-0.5);

	int i = ix + nax_*(iy + nay_*iz);
	
	double r2 = x2 + y2 + z2;

	values_[i] = (r2 < R2) ? 1.0 : 0.0;
      }
    }
  }

}

//----------------------------------------------------------------------

void CharmPatch::p_receive(int axis, int face, int n, double * buffer_ghost)
{
  // Update ghosts on given face from buffer

  buffer_to_ghost_(axis, face, buffer_ghost);

  
  if (receives_.remaining() == 0) {
    compute_();
  }
}

//----------------------------------------------------------------------
void CharmPatch::compute_()
{
  const int dx = 1;
  const int dy = nax_;
  const int dz = nax_*nay_;

  const double o6 = 1.0 / 6.0;

  double * values_old = new double [nax_*nay_*naz_];
  for (int i=0; i<nax_*nay_*naz_; i++) values_old[i] = values_[i];

  for (int iz=izl_; iz<izu_; iz++) {
    for (int iy=iyl_; iy<iyu_; iy++) {
      for (int ix=ixl_; ix<ixu_; ix++) {
	int i = ix + nax_*(iy + nay_*iz);
	values_[i] = o6*(values_old[i-dx] + values_old[i+dx] +
			 values_old[i-dy] + values_old[i+dy] +
			 values_old[i-dz] + values_old[i+dz]);
      }
    }
  }

  delete [] values_old;

  write_();

  main_proxy_.p_next(thisIndex.x,thisIndex.y,thisIndex.z,norm_());

  // 
}

//----------------------------------------------------------------------

void CharmPatch::p_evolve() 
{

  // Update ghost zones

  face_to_buffer_(0,0,buffer_[0][0]);
  face_to_buffer_(0,1,buffer_[0][1]);
  face_to_buffer_(1,0,buffer_[1][0]);
  face_to_buffer_(1,1,buffer_[1][1]);
  face_to_buffer_(2,0,buffer_[2][0]);
  face_to_buffer_(2,1,buffer_[2][1]);
      
  receive_();

  // (...execution continued in last receive executed locally)
}

//======================================================================

void CharmPatch::receive_()
{
  CProxy_CharmPatch patch = thisProxy;

  int ix = thisIndex.x;
  int iy = thisIndex.y;
  int iz = thisIndex.z;

  int ixm = (ix - 1 + nbx_) % nbx_;
  int iym = (iy - 1 + nby_) % nby_;
  int izm = (iz - 1 + nbz_) % nbz_;
  int ixp = (ix + 1) % nbx_;
  int iyp = (iy + 1) % nby_;
  int izp = (iz + 1) % nbz_;

  int nx = ngx_*nay_*naz_;
  int ny = nax_*ngy_*naz_;
  int nz = nax_*nay_*ngz_;

  patch(ixp,iy,iz).p_receive(0,0,nx,buffer_[0][1]);
  patch(ixm,iy,iz).p_receive(0,1,nx,buffer_[0][0]);
  patch(ix,iyp,iz).p_receive(1,0,ny,buffer_[1][1]);
  patch(ix,iym,iz).p_receive(1,1,ny,buffer_[1][0]);
  patch(ix,iy,izp).p_receive(2,0,nz,buffer_[2][1]);
  patch(ix,iy,izm).p_receive(2,1,nz,buffer_[2][0]);
}

//======================================================================

bool CharmPatch::clear_boundary_(int axis, int face, double * buffer)
{
  if ((axis == 0 && face == 0 && thisIndex.x == 0) ||
      (axis == 0 && face == 1 && thisIndex.x == nbx_-1)) {
    for (int i=0; i<ngx_*nay_*naz_; i++) {
      buffer[i] = 0.0;
    }
    return true;
  } else if ((axis == 1 && face == 0 && thisIndex.y == 0) ||
	     (axis == 1 && face == 1 && thisIndex.y == nby_-1)) {
    for (int i=0; i<nax_*ngy_*naz_; i++) {
      buffer[i] = 0.0;
    }
    return true;
  } else if ((axis == 2 && face == 0 && thisIndex.z== 0) ||
	     (axis == 2 && face == 1 && thisIndex.z == nbz_-1)) {
    for (int i=0; i<nax_*nay_*ngz_; i++) {
      buffer[i] = 0.0;
    }
    return true;
  } else {
    return false;
  }
}

//======================================================================

void CharmPatch::face_to_buffer_(int axis, int face, double * buffer)
{

  if (clear_boundary_(axis,face,buffer)) return;

  int n  = nax_*nay_*naz_;

  int nx = ngx_*nay_*naz_;
  int ny = nax_*ngy_*naz_;
  int nz = nax_*nay_*ngz_;

  int kdx = face*(nvx_-ngx_) + ngx_;
  int kdy = face*(nvy_-ngy_) + ngy_;
  int kdz = face*(nvz_-ngz_) + ngz_;

  if (axis == 0) {
    for (int iz=0; iz<naz_; iz++) {
      for (int iy=0; iy<nay_; iy++) {
	for (int ix=0; ix<ngx_; ix++) {
	  int kx = kdx +ix;
	  int iv = kx + nax_*(iy + nay_*iz);
	  int ig = iy + nay_*(iz + naz_*ix);
	  assert(0 <= ig && ig < nx);
	  assert(0 <= iv && iv < n);
	  buffer[ig] = values_[iv];
	}
      }
    }
  } else if (axis == 1) {
    for (int iz=0; iz<naz_; iz++) {
      for (int iy=0; iy<ngy_; iy++) {
	for (int ix=0; ix<nax_; ix++) {
	  int ky = kdy + iy;
	  int iv = ix + nax_*(ky + nay_*iz);
	  int ig = iz + naz_*(ix + nax_*iy);
	  assert(0 <= ig && ig < ny);
	  assert(0 <= iv && iv < n);
	  buffer[ig] = values_[iv];
	}
      }
    }
  } else if (axis == 2) {
    for (int iz=0; iz<ngz_; iz++) {
      for (int iy=0; iy<nay_; iy++) {
	for (int ix=0; ix<nax_; ix++) {
	  int kz = kdz + iz;
	  int iv = ix + nax_*(iy + nay_*kz);
	  int ig = ix + nax_*(iy + nay_*iz);
	  assert(0 <= ig && ig < nz);
	  assert(0 <= iv && iv < n);
	  buffer[ig] = values_[iv];
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void CharmPatch::buffer_to_ghost_(int axis, int face, double * buffer)
{
  clear_boundary_(axis,face,buffer);

  int n  = nax_*nay_*naz_;

  int nx = ngx_*nay_*naz_;
  int ny = nax_*ngy_*naz_;
  int nz = nax_*nay_*ngz_;

  int kdx = face*(nax_-ngx_);
  int kdy = face*(nay_-ngy_);
  int kdz = face*(naz_-ngz_);

  if (axis == 0) {
    for (int iz=0; iz<naz_; iz++) {
      for (int iy=0; iy<nay_; iy++) {
	for (int ix=0; ix<ngx_; ix++) {
	  int kx = kdx + ix;
	  int iv = kx + nax_*(iy + nay_*iz);
	  int ig = iy + nay_*(iz + naz_*ix);
	  assert(0 <= ig && ig < nx);
	  assert(0 <= iv && iv < n);
	  values_[iv] = buffer[ig];
	}
      }
    }
  } else if (axis == 1) {
    for (int iz=0; iz<naz_; iz++) {
      for (int iy=0; iy<ngy_; iy++) {
	for (int ix=0; ix<nax_; ix++) {
	  int ky = kdy + iy;
	  int iv = ix + nax_*(ky + nay_*iz);
	  int ig = iz + naz_*(ix + nax_*iy);
	  assert(0 <= ig && ig < ny);
	  assert(0 <= iv && iv < n);
	  values_[iv] = buffer[ig];
	}
      }
    }
  } else if (axis == 2) {
    for (int iz=0; iz<ngz_; iz++) {
      for (int iy=0; iy<nay_; iy++) {
	for (int ix=0; ix<nax_; ix++) {
	  int kz = kdz + iz;
	  int iv = ix + nax_*(iy + nay_*kz);
	  int ig = ix + nax_*(iy + nay_*iz);
	  assert(0 <= ig && ig < nz);
	  assert(0 <= iv && iv < n);
	  values_[iv] = buffer[ig];
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void CharmPatch::print_ ()
{
  int ix0 = thisIndex.x * nvx_;
  int iy0 = thisIndex.y * nvy_;
  int iz0 = thisIndex.z * nvz_;
  for (int iz=izl_; iz<izu_; iz++) {
    for (int iy=iyl_; iy<iyu_; iy++) {
      for (int ix=ixl_; ix<ixu_; ix++) {
	int i = ix + nax_*(iy + nay_*iz);
	CkPrintf ("%03d %03d %03d %g\n",ix+ix0,iy+iy0,iz+iz0,values_[i]);
      }
    }
  }
}

//----------------------------------------------------------------------

double CharmPatch::norm_ ()
{
  double s2 = 0.0;
  for (int iz=izl_; iz<izu_; iz++) {
    for (int iy=iyl_; iy<iyu_; iy++) {
      for (int ix=ixl_; ix<ixu_; ix++) {
	int i = ix + nax_*(iy + nay_*iz);
	s2 += values_[i]*values_[i];
      }
    }
  }
  return s2;
}

//----------------------------------------------------------------------

void CharmPatch::write_ ()
{
  if (cycle_values_ % cycle_write_ == 0) {

    hid_t file;
    hid_t dataset;
    hid_t dataspace;
    hid_t datatype;
    const hsize_t nv[3] = {nax_,nay_,naz_};
    const hsize_t na[3] = {nax_,nay_,naz_};
    //    int i0 = ngx_ + nax_*(ngy_ + nay_*ngz_);
    int i0 = 0;

    datatype = H5T_NATIVE_DOUBLE;
    char filename[80];

    sprintf(filename,"jacobi-%05d-%02d-%02d-%02d.hdf5",
	    cycle_values_,thisIndex.x,thisIndex.y,thisIndex.z);
    file      = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert (file > 0);
    dataspace = H5Screate_simple (3,nv,na);
    assert (dataspace > 0);
    dataset   = H5Dcreate( file, "dataset",datatype,dataspace, 
			   H5P_DEFAULT );
    assert (dataset > 0);
    H5Dwrite (dataset, datatype, dataspace, H5S_ALL, H5P_DEFAULT, values_+i0);
    H5Dclose( dataset);
    H5Sclose (dataspace);
    H5Fclose (file);
  }
}
