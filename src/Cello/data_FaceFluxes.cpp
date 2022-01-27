// See LICENSE_CELLO file for license and copyright information

/// @file     data_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    [\ref Data] Implementation of the FaceFluxes class

#include "data.hpp"

// #define TRACE_FACE_FLUXES

#ifdef TRACE_FACE_FLUXES
#   undef TRACE_FACE_FLUXES
#   define TRACE_FACE_FLUXES(MSG) CkPrintf ("%d TRACE_FACE_FLUXES %s\n",CkMyPe(),MSG);
#else
#   define TRACE_FACE_FLUXES(MSG) /* ... */
#endif

//----------------------------------------------------------------------

FaceFluxes::FaceFluxes (Face face, int index_field,
                        int nx, int ny, int nz,
                        int cx, int cy, int cz)
  : face_(face),
    index_field_(index_field),
    nx_(nx),ny_(ny),nz_(nz),
    cx_(cx),cy_(cy),cz_(cz),
    delete_fluxes_(false),
    fluxes_(nullptr)
{
  int fx,fy,fz; // face
  face.get_face(&fx,&fy,&fz);
  // set face size to 1 when orthogonal to face
  // (should be 0 for axes greater than rank)
  if (fx && nx) nx_ = 1;
  if (fy && ny) ny_ = 1;
  if (fz && nz) nz_ = 1;
  ASSERT3 ("FaceFluxes::FaceFluxes",
          "Block dimensions %d %d %d must be even or 1",
          nx_,ny_,nz_,
           (( nx_==1 || ((nx_&1)==0) ) &&
            ( ny_==1 || ((ny_&1)==0) ) &&
            ( nz_==1 || ((nz_&1)==0) )));
}

//----------------------------------------------------------------------

void FaceFluxes::pup (PUP::er &p)
{
  TRACEPUP;

  p | face_;

  p | index_field_;

  p | nx_;
  p | ny_;
  p | nz_;

  p | cx_;
  p | cy_;
  p | cz_;

  p | delete_fluxes_;
  
  int fluxes_allocated = (fluxes_ != nullptr);
  p | fluxes_allocated;
  if (fluxes_allocated) {
    if (delete_fluxes_) {
      // depends on nx,ny,nz,cx,cy,cz
      int size = get_size();
      p | size;
      if (p.isUnpacking()) allocate_storage();
      PUParray(p,fluxes_,size);
    } else {
      // all Block fluxes_ points into single array: PUP values
      // handled by parent FluxData object
    }
  } else {
    fluxes_ = nullptr;
  }
}

//----------------------------------------------------------------------

void FaceFluxes::set_flux_array ( std::vector<cello_float> array,
                                  int dx, int dy, int dz)
{
  TRACE_FACE_FLUXES("set_flux_array()");
  int mx,my,mz;
  get_size(&mx,&my,&mz);
  ASSERT5("FaceFluxes::set_flux_array",
          "Input array size %lu is smaller than required %d = %d * %d *%d",
          array.size(),mx*my*mz,mx,my,mz,
          ((int)array.size() >= mx*my*mz) );
  ASSERT4("FaceFluxes::set_flux_array",
          "Input array size %lu is too small for strides (dx,dy,dz)=(%d,%d,%d)",
          array.size(),dx,dy,dz,
          ((int)array.size() >= (mx-1)*dx+(my-1)*dy+(mz-1)*dz) );

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        int i_f = ix + mx*(iy+my*iz);
        int i_a = ix*dx + iy*dy + iz*dz;
        fluxes_[i_f] = array[i_a];
      }
    }
  }
}
  
//----------------------------------------------------------------------

cello_float * FaceFluxes::flux_array (int * dx, int * dy, int *dz)
{
  TRACE_FACE_FLUXES("flux_array()");
  int mx,my,mz;
  get_size(&mx,&my,&mz);
  if (dx) (*dx) = 1;
  if (dy) (*dy) = mx;
  if (dz) (*dz) = mx*my;

  return fluxes_;
}
  
//----------------------------------------------------------------------

void FaceFluxes::coarsen (int cx, int cy, int cz, int rank)
{
  TRACE_FACE_FLUXES("coarsen()");
  int mxf,myf,mzf;
  int m = get_size(&mxf,&myf,&mzf);
  
  cello_float * fluxes_fine = new cello_float[m];
  std::copy_n(fluxes_,m,fluxes_fine);
  std::fill_n(fluxes_,m,0.0);

  const int mxc = (nx_ > 1) ? (nx_/2 + cx_) : 1;
  const int myc = (ny_ > 1) ? (ny_/2 + cy_) : 1;
  const int mzc = (nz_ > 1) ? (nz_/2 + cz_) : 1;

  if (rank == 2) {

    const int dxc = 1;
    const int dyc = nx_;
    
    const int ixc0 = cx*mxc*dxc;
    const int iyc0 = cy*myc*dyc;
    
    const cello_float dArea = 0.25;
    if (mxf == 1) {

      const int ixf = 0;
      const int ic0 = iyc0;
      for (int iyf=0; iyf<myf; iyf++) {
        int i_f = ixf + mxf*iyf;
        int iycm = dyc*(iyf>>1);
        int iycp = dyc*((iyf+cy_)>>1);
        fluxes_[ic0+iycm] += 0.5*dArea*fluxes_fine[i_f];
        fluxes_[ic0+iycp] += 0.5*dArea*fluxes_fine[i_f];
      }

    } else if (myf == 1) {

      const int iyf = 0;
      const int ic0 = ixc0;
      for (int ixf=0; ixf<mxf; ixf++) {
        int i_f = ixf + mxf*iyf;
        int ixcm = dxc*(ixf>>1);
        int ixcp = dxc*((ixf+cx_)>>1);
        fluxes_[ic0+ixcm] += 0.5*dArea*fluxes_fine[i_f];
        fluxes_[ic0+ixcp] += 0.5*dArea*fluxes_fine[i_f];
      }
    }

  } else if (rank == 3) {
    
    const int dxc = 1;
    const int dyc = nx_;
    const int dzc = nx_*ny_;
    
    const int ixc0 = cx*mxc*dxc;
    const int iyc0 = cy*myc*dyc;
    const int izc0 = cz*mzc*dzc;

    const cello_float dArea = 0.125;
    
    if (mxf == 1) {
      
      const int ixf = 0;
      const int ic0 =iyc0 + izc0;
      for (int izf=0; izf<mzf; izf++) {
        int izcm = dzc*(izf>>1);
        int izcp = dzc*((izf+cz_)>>1);
        for (int iyf=0; iyf<myf; iyf++) {
          int iycm = dyc*(iyf>>1);
          int iycp = dyc*((iyf+cy_)>>1);
          int i_f = ixf + mxf*(iyf + myf*izf);
          fluxes_[ic0+izcm+iycm] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+izcm+iycp] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+izcp+iycm] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+izcp+iycp] += 0.25*dArea*fluxes_fine[i_f];
        }
      }
    } else if (myf == 1) {

      const int iyf=0;
      const int ic0 =ixc0 + izc0;
      for (int izf=0; izf<mzf; izf++) {
        int izcm = dzc*(izf>>1);
        int izcp = dzc*((izf+cz_)>>1);
        for (int ixf=0; ixf<mxf; ixf++) {
          int ixcm = dxc*(ixf>>1);
          int ixcp = dxc*((ixf+cx_)>>1);
          int i_f = ixf + mxf*(iyf + myf*izf);
          fluxes_[ic0+izcm+ixcm] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+izcm+ixcp] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+izcp+ixcm] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+izcp+ixcp] += 0.25*dArea*fluxes_fine[i_f];
        }
      }
    } else if (mzf == 1) {

      const int izf=0;
      const int ic0 =ixc0 + iyc0;
      for (int iyf=0; iyf<myf; iyf++) {
        int iycm = dyc*(iyf>>1);
        int iycp = dyc*((iyf+cy_)>>1);
        for (int ixf=0; ixf<mxf; ixf++) {
          int ixcm = dxc*(ixf>>1);
          int ixcp = dxc*((ixf+cx_)>>1);
          int i_f = ixf + mxf*(iyf + myf*izf);
          fluxes_[ic0+iycm+ixcm] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+iycm+ixcp] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+iycp+ixcm] += 0.25*dArea*fluxes_fine[i_f];
          fluxes_[ic0+iycp+ixcp] += 0.25*dArea*fluxes_fine[i_f];
        }
      }
    }

  } else {

    ERROR1 ("FaceFluxes::coarsen()",
            "Unsupported rank = %d",
            rank);

  }

  delete [] fluxes_fine;
}
  
//----------------------------------------------------------------------

void FaceFluxes::accumulate
(const FaceFluxes & ff_2, int cx,int cy, int cz, int rank)
{
  TRACE_FACE_FLUXES("accumulate()");
  FaceFluxes & ff_1 = *this;

  int mx1,my1,mz1;
  int mx2,my2,mz2;
  
  ff_1.get_size(&mx1,&my1,&mz1);
  ff_2.get_size(&mx2,&my2,&mz2);

  for (int iz=0; iz<mz1; iz++) {
    for (int iy=0; iy<my1; iy++) {
      for (int ix=0; ix<mx1; ix++) {
        int i1=ix + mx1*(iy + my1*iz);
        int i2=ix + mx2*(iy + my2*iz);
        ff_1.fluxes_[i1] += ff_2.fluxes_[i2];
      }
    }
  }
 
}
  
//----------------------------------------------------------------------

FaceFluxes & FaceFluxes::operator *= (double weight)
{
  TRACE_FACE_FLUXES("operator *=()");
  int mx,my,mz;
  get_size(&mx,&my,&mz);
  const int m=mx*my*mz;

  for (int i=0; i<m; i++) fluxes_[i] *= weight;
  
  return *this;
}
  
//======================================================================

int FaceFluxes::data_size () const
{
  TRACE_FACE_FLUXES("data_size()");
  int size = 0;

  size += face_.data_size();
  SIZE_SCALAR_TYPE(size,int,index_field_);
  SIZE_SCALAR_TYPE(size,int,nx_);
  SIZE_SCALAR_TYPE(size,int,ny_);
  SIZE_SCALAR_TYPE(size,int,nz_);
  SIZE_SCALAR_TYPE(size,int,cx_);
  SIZE_SCALAR_TYPE(size,int,cy_);
  SIZE_SCALAR_TYPE(size,int,cz_);
  SIZE_SCALAR_TYPE(size,int,delete_fluxes_);
  
  SIZE_ARRAY_TYPE(size,cello_float,fluxes_,get_size());

  return size;
}

//----------------------------------------------------------------------

char * FaceFluxes::save_data (char * buffer) const
{
  TRACE_FACE_FLUXES("save_data()");
  union {
    int  * pi;
    char * pc;
    double * pd;
    cello_float * pcf;
  };

  pc = (char *) buffer;

  pc = face_.save_data(pc);

  SAVE_SCALAR_TYPE(pc,int,index_field_);
  SAVE_SCALAR_TYPE(pc,int,nx_);
  SAVE_SCALAR_TYPE(pc,int,ny_);
  SAVE_SCALAR_TYPE(pc,int,nz_);
  SAVE_SCALAR_TYPE(pc,int,cx_);
  SAVE_SCALAR_TYPE(pc,int,cy_);
  SAVE_SCALAR_TYPE(pc,int,cz_);
  SAVE_SCALAR_TYPE(pc,int,delete_fluxes_);

  SAVE_ARRAY_TYPE(pc,cello_float,fluxes_,get_size());
  
  ASSERT2("FaceFluxes::save_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------

char * FaceFluxes::load_data (char * buffer)
{
  TRACE_FACE_FLUXES("load_data()");
  union {
    int  * pi;
    char * pc;
    double * pd;
    cello_float * pcf;
  };

  pc = (char *) buffer;

  pc = face_.load_data(pc);
  
  LOAD_SCALAR_TYPE(pc,int,index_field_);
  LOAD_SCALAR_TYPE(pc,int,nx_);
  LOAD_SCALAR_TYPE(pc,int,ny_);
  LOAD_SCALAR_TYPE(pc,int,nz_);
  LOAD_SCALAR_TYPE(pc,int,cx_);
  LOAD_SCALAR_TYPE(pc,int,cy_);
  LOAD_SCALAR_TYPE(pc,int,cz_);
  LOAD_SCALAR_TYPE(pc,int,delete_fluxes_);
  delete_fluxes_ = true;
  fluxes_ = new cello_float [get_size()];
  LOAD_ARRAY_TYPE(pc,cello_float,fluxes_,get_size());
  ASSERT2("FaceFluxes::load_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------

void FaceFluxes::print (Block * block, std::string message)
{
  cello_float min=1e30, max=-1e30, avg=0.0;
  int mx,my,mz;
  const int n = get_size(&mx,&my,&mz);
  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        const int i = ix + mx*(iy+my*iz);
        min = std::min(min,fluxes_[i]);
        avg += fluxes_[i];
        max = std::max(max,fluxes_[i]);
      }
    }
  }
  avg /= mx*my*mz;
  CkPrintf ("FaceFluxes::print %s %p %s %d %d %d %g %g %g\n",
            block->name().c_str(),(void *)fluxes_,
            message.c_str(),mx,my,mz,min,avg,max);
}
