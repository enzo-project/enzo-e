// See LICENSE_CELLO file for license and copyright information

/// @file     data_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    [\ref Data] Implementation of the FaceFluxes class

#include "data.hpp"

// #define DEBUG_FACE_FLUXES
// #define DEBUG_REFRESH

FaceFluxes::FaceFluxes (Face face, int index_field,
                        int nx, int ny, int nz,
                        int level, double dt,
                        int cx, int cy, int cz)
  : face_(face),
    level_block_(level),
    level_neighbor_(level),
    dt_block_(dt),
    dt_neighbor_(dt),
    fluxes_(),
    index_field_(index_field),
    nx_(nx),ny_(ny),nz_(nz),
    cx_(cx),cy_(cy),cz_(cz)
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

  p | level_block_;
  p | level_neighbor_;
  p | dt_block_;
  p | dt_neighbor_;

  p | fluxes_;

  p | index_field_;

  p | nx_;
  p | ny_;
  p | nz_;

  p | cx_;
  p | cy_;
  p | cz_;
}

//----------------------------------------------------------------------

void FaceFluxes::set_flux_array ( std::vector<double> array,
                                  int dx, int dy, int dz)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  ASSERT5("FaceFluxes::set_flux_array",
          "Input array size %lu is smaller than required %d = %d * %d *%d",
          array.size(),mx*my*mz,mx,my,mz,
          (array.size() >= mx*my*mz) );
  ASSERT4("FaceFluxes::set_flux_array",
          "Input array size %lu is too small for strides (dx,dy,dz)=(%d,%d,%d)",
          array.size(),dx,dy,dz,
          (array.size() >= (mx-1)*dx+(my-1)*dy+(mz-1)*dz) );

  ASSERT5("FaceFluxes::set_flux_array",
          "Flux array size %lu is smaller than required %d = %d * %d *%d",
          fluxes_.size(),mx*my*mz,mx,my,mz,
          (fluxes_.size() >= mx*my*mz) );

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

std::vector<double> & FaceFluxes::flux_array (int * dx, int * dy, int *dz)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  if (dx) (*dx) = 1;
  if (dy) (*dy) = mx;
  if (dz) (*dz) = mx*my;

  return fluxes_;
}
  
//----------------------------------------------------------------------

float ratio_cell_volume
(const FaceFluxes & ff_1, const FaceFluxes & ff_2, const int rank)
{
  const int L_1 = ff_1.level_neighbor();
  const int L_2 = ff_2.level_neighbor();
  if (L_1 <= L_2) {
    const int shift = (1 << (L_2 - L_1)*rank);
    return 1.0 * shift;
  } else if (L_1 > L_2) {
    const int shift = (1 << (L_1 - L_2)*rank);
    return 1.0 / shift;
  }
}

//----------------------------------------------------------------------

float ratio_time_step (const FaceFluxes & ff_1, const FaceFluxes & ff_2)
{
  return ff_1.time_step_block() / ff_2.time_step_block();
}
  
//----------------------------------------------------------------------

void FaceFluxes::coarsen ()
{
  std::vector<double> fluxes_fine = fluxes_;
  int mxf,myf,mzf;
  get_dimensions(&mxf,&myf,&mzf);
  // not correct??  need nxf/2 for ghost and centering

  // nxc + a = nx/2 + a
  const int mxc = (nx_ > 1) ? (nx_/2 + cx_) : 1;
  const int myc = (ny_ > 1) ? (ny_/2 + cy_) : 1;
  const int mzc = (nz_ > 1) ? (nz_/2 + cz_) : 1;

  fluxes_.resize(mxc*myc*mzc);

  std::fill(fluxes_.begin(),fluxes_.end(),0.0);

  const int dxc = 1;
  const int dyc = mxc;
  const int dzc = mxc*myc;
  if (mxf == 1) {
    int ixf=0;
    for (int izf=0; izf<mzf; izf++) {
      int izcm = izf>>1;
      int izcp = (izf+cz_)>>1;
      for (int iyf=0; iyf<myf; iyf++) {
        int iycm = iyf>>1;
        int iycp = (iyf+cy_)>>1;
        int i_f = ixf + mxf*(iyf + myf*izf);
        fluxes_[iycm*dyc + izcm*dzc] += fluxes_fine[i_f];
        fluxes_[iycm*dyc + izcp*dzc] += fluxes_fine[i_f];
        fluxes_[iycp*dyc + izcm*dzc] += fluxes_fine[i_f];
        fluxes_[iycp*dyc + izcp*dzc] += fluxes_fine[i_f];
      }
    }
  } else if (myf == 1) {
    int iyf=0;
    for (int izf=0; izf<mzf; izf++) {
      int izcm = izf>>1;
      int izcp = (izf+cz_)>>1;
      for (int ixf=0; ixf<mxf; ixf++) {
        int ixcm = ixf>>1;
        int ixcp = (ixf+cx_)>>1;
        int i_f = ixf + mxf*(iyf + myf*izf);
        fluxes_[izcm*dzc + ixcm*dxc] += fluxes_fine[i_f];
        fluxes_[izcp*dzc + ixcm*dxc] += fluxes_fine[i_f];
        fluxes_[izcm*dzc + ixcp*dxc] += fluxes_fine[i_f];
        fluxes_[izcp*dzc + ixcp*dxc] += fluxes_fine[i_f];
      }
    }
  } else if (mzf == 1) {
    int izf=0;
    for (int iyf=0; iyf<myf; iyf++) {
      int iycm = iyf>>1;
      int iycp = (iyf+cy_)>>1;
      for (int ixf=0; ixf<mxf; ixf++) {
        int ixcm = ixf>>1;
        int ixcp = (ixf+cx_)>>1;
        int i_f = ixf + mxf*(iyf + myf*izf);
        fluxes_[ixcm*dxc + iycm*dyc] += fluxes_fine[i_f];
        fluxes_[ixcp*dxc + iycm*dyc] += fluxes_fine[i_f];
        fluxes_[ixcm*dxc + iycp*dyc] += fluxes_fine[i_f];
        fluxes_[ixcp*dxc + iycp*dyc] += fluxes_fine[i_f];
      }
    }
  }
  nx_ = (nx_ != 1) ? nx_>>1 : 1;
  ny_ = (ny_ != 1) ? ny_>>1 : 1;
  nz_ = (nz_ != 1) ? nz_>>1 : 1;
  --level_neighbor_;
}
  
//----------------------------------------------------------------------

void FaceFluxes::accumulate (const FaceFluxes & ff_2, int cx,int cy, int cz)
{
  FaceFluxes & ff_1 = *this;

  int ix01,iy01,iz01;
  int nx1,ny1,nz1;
  int mx1,my1,mz1;
  
  ff_1.get_start(&ix01,&iy01,&iz01,cx,cy,cz,ff_2);
  ff_1.get_size(&nx1,&ny1,&nz1,ff_2);
  ff_1.get_dimensions(&mx1,&my1,&mz1);
  
  int ix02,iy02,iz02;
  int nx2,ny2,nz2;
  int mx2,my2,mz2;
  ff_2.get_start(&ix02,&iy02,&iz02,cx,cy,cz,ff_1);
  ff_2.get_size(&nx2,&ny2,&nz2,ff_1);
  ff_2.get_dimensions(&mx2,&my2,&mz2);

  ASSERT6 ("FaceFluxes::operator +=()",
          "Flux array sizes (%d %d %d) and (%d %d %d) must be conforming",
           nx1,ny1,nz1,nx2,ny2,nz2,
           ((nx1==nx2) && (ny1==ny2) && (nz1==nz2)));

  for (int iz=0; iz<nz1; iz++) {
    for (int iy=0; iy<ny1; iy++) {
      for (int ix=0; ix<nx1; ix++) {
        int i1=(ix+ix01) + mx1*((iy+iy01) + my1*(iz+iz01));
        int i2=(ix+ix02) + mx2*((iy+iy02) + my2*(iz+iz02));
        ff_1.fluxes_[i1] += ff_2.fluxes_[i2];
      }
    }
  }
}
  
//----------------------------------------------------------------------

FaceFluxes & FaceFluxes::operator *= (double weight)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  const int m=mx*my*mz;

  for (int i=0; i<m; i++) fluxes_[i] *= weight;
  
  return *this;
}
  
//======================================================================

int FaceFluxes::data_size () const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::data_size()\n",CkMyPe());
#endif  
  int size = 0;

  size += face_.data_size();
  size += sizeof(int);    // std::vector<double> fluxes_;
  int n = fluxes_.size();
  size += n*sizeof(double);
  size += sizeof(int);    // int level_block_;
  size += sizeof(int);    // int level_neighbor_;
  size += sizeof(double);    // double dt_block_;
  size += sizeof(double);    // double dt_neighbor_;
  size += sizeof(int);    // int index_field_;
  size += 3*sizeof(int);  // int nx_,ny_,nz_;
  size += 3*sizeof(int);  // int cx_,cy_,cz_;

  return size;
}

//----------------------------------------------------------------------

char * FaceFluxes::save_data (char * buffer) const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::save_data()\n",CkMyPe());
#endif  
  union {
    int  * pi;
    char * pc;
    double * pd;
  };

  pc = (char *) buffer;

  pc = face_.save_data(pc);
  int n = fluxes_.size();
  (*pi++) = n;
  for (int i=0; i<n; i++) {
    (*pd++) = fluxes_[i];
  }
  (*pi++) = level_block_;
  (*pi++) = level_neighbor_;
  (*pd++) = dt_block_;
  (*pd++) = dt_neighbor_;
  (*pi++) = index_field_;
  (*pi++) = nx_;
  (*pi++) = ny_;
  (*pi++) = nz_;
  (*pi++) = cx_;
  (*pi++) = cy_;
  (*pi++) = cz_;

  ASSERT2("FaceFluxes::save_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

//----------------------------------------------------------------------

char * FaceFluxes::load_data (char * buffer)
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::load_data()\n",CkMyPe());
#endif  
  union {
    int  * pi;
    char * pc;
    double * pd;
  };

  pc = (char *) buffer;

  pc = face_.load_data(pc);
  int n = (*pi++);
  fluxes_.resize(n);
  for (int i=0; i<n; i++) {
    fluxes_[i] = (*pd++);
  }
  level_block_ = (*pi++);
  level_neighbor_ = (*pi++);
  dt_block_ = (*pd++);
  dt_neighbor_ = (*pd++);
  index_field_ = (*pi++);
  nx_ = (*pi++);
  ny_ = (*pi++);
  nz_ = (*pi++);
  cx_ = (*pi++);
  cy_ = (*pi++);
  cz_ = (*pi++);

  ASSERT2("FaceFluxes::load_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

  return pc;
}

