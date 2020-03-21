// See LICENSE_CELLO file for license and copyright information

/// @file     data_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    [\ref Data] Implementation of the FaceFluxes class

#include "data.hpp"

// #define DEBUG_FACE_FLUXES

FaceFluxes::FaceFluxes (Face face, int index_field,
                        int nx, int ny, int nz,
                        double hx, double hy, double hz, double dt)
  : face_(face),
    index_field_(index_field),
    fluxes_(),
    hx_(hx),hy_(hy),hz_(hz),
    dt_(dt),
    cx_(0),cy_(0),cz_(0),
    gx_(0),gy_(0),gz_(0)
{
  bool lx,ly,lz;
  face_.adjacency(&lx,&ly,&lz);
  nx_ = lx ? nx : 1;
  ny_ = ly ? ny : 1;
  nz_ = lz ? nz : 1;
  update_dimensions_();
  ASSERT3 ("FaceFluxes::FaceFluxes",
          "Block dimensions must be even %d %d %d",
          nx,ny,nz,
          ( ((nx&1)==0) &&   // even number of x-axis cells
            (((ny&1)==0) || (ny==1 && nz==1)) && // even y-axis or rank==1
            (((nz&1)==0) || (nz==1))) );         // even z-axis or rank==2
}

//----------------------------------------------------------------------

void FaceFluxes::pup (PUP::er &p)
{
  TRACEPUP;

  p | face_;
  p | index_field_;

  p | mx_;
  p | my_;
  p | mz_;

  p | nx_;
  p | ny_;
  p | nz_;

  p | cx_;
  p | cy_;
  p | cz_;

  p | gx_;
  p | gy_;
  p | gz_;

  p | fluxes_;
    
  p | hx_;
  p | hy_;
  p | hz_;
    
  p | dt_;

}

//----------------------------------------------------------------------

void FaceFluxes::set_ghost(int gx, int gy, int gz)
{
  gx_ = gx;
  gy_ = gy;
  gz_ = gz;
  update_dimensions_();
}

//----------------------------------------------------------------------

void FaceFluxes::set_centering(int cx, int cy, int cz)
{
  cx_ = cx;
  cy_ = cy;
  cz_ = cz;
  update_dimensions_();
}
  
//----------------------------------------------------------------------

void FaceFluxes::allocate ()
{
  int mx,my,mz;
  get_dimensions (&mx,&my,&mz);

  const int rank = face_.rank();
  const int m = (rank == 1) ? mx : (rank == 2) ? mx*my : mx*my*mz;
  
  ASSERT4 ("FaceFluxes::allocate()",
          "Unexpected flux dimensions %d %d %d for rank %d",
           mx,my,mz,rank, (m > 0) );
  
  fluxes_.resize(m);
}

//----------------------------------------------------------------------

void FaceFluxes::deallocate()
{
  fluxes_.resize(0);
}

//----------------------------------------------------------------------

void FaceFluxes::clear()
{
  std::fill (fluxes_.begin(),fluxes_.end(),0.0);
}

//----------------------------------------------------------------------

Face FaceFluxes::face () const
{
  return face_;
}

//----------------------------------------------------------------------

void FaceFluxes::get_element_size (double *hx, double *hy, double * hz) const
{
  if (hx) (*hx) = hx_;
  if (hy) (*hy) = hy_;
  if (hz) (*hz) = hz_;
}
  
//----------------------------------------------------------------------

void FaceFluxes::get_limits
(int * ixl, int * ixu,
 int * iyl, int * iyu,
 int * izl, int * izu) const
{
  /// If same level, will be
  ///
  /// (ixl,iyl,izl) == (0,0,0) and
  /// (ixl,iyl,izl) == (mx,my,mz);
  ///
  /// if neighbor is finer will depend on child bits:
  ///
  /// (ixl,iyl,izl) == (cx*mx/2,cy*my/2,cz*mz/2)
  /// (ixu,iyu,izu) == ((2-cx)*mx/2,(2-cy)*my/2,(2-cz)*mz/2)
  ///
  /// o-----------o
  /// | *-------* |
  /// | |       | |
  /// | |       *---*  iyu
  /// | |       |   |
  /// | *-------*---*  iyl
  /// o-----------o    0


  int level_block    = face().index_block().level();
  int level_neighbor = face().index_neighbor().level();
  
  // ASSERT2("FaceFluxes::get_limits",
  //         "Assumes neighbor level %d is no coarser than block level %d",
  //         level_block,level_neighbor,
  //         (level_block <= level_neighbor));

#ifdef DEBUG_FACE_FLUXES
  CkPrintf ("DEBUG_FLUX level (block neighbor) = (%d %d)\n",
            level_block,level_neighbor);
#endif
  
  if (level_block == level_neighbor) {
    if (ixl != nullptr) (*ixl) = gx_;
    if (ixu != nullptr) (*ixu) = nx_+gx_;
    if (iyl != nullptr) (*iyl) = gy_;
    if (iyu != nullptr) (*iyu) = ny_+gy_;
    if (izl != nullptr) (*izl) = gz_;
    if (izu != nullptr) (*izu) = nz_+gz_;
  } else if (level_block - level_neighbor == 1) {
    // neighbor is coarser
    INCOMPLETE ("FaceFluxes::get_limits() coarse neighbor");
  } else if (level_neighbor - level_block == - 1) {
    // neighbor is finer
    INCOMPLETE ("FaceFluxes::get_limits() fine neighbor");
  }
    
  /// (ixl,iyl,izl) == (0,0,0) and
  /// (ixl,iyl,izl) == (mx,my,mz);
  /// if neighbor is finer will depend on child bits:
  /// (ixl,iyl,izl) == (cx*mx/2,cy*my/2,cz*mz/2)
  /// (ixu,iyu,izu) == ((2-cx)*mx/2,(2-cy)*my/2,(2-cz)*mz/2)
}

//----------------------------------------------------------------------

double FaceFluxes::time_step () const
{
  return dt_;
}
  
//----------------------------------------------------------------------

void FaceFluxes::get_dimensions (int *mx, int *my, int *mz) const
{
  const int rank = face_.rank();

  if (mx) (*mx) = (rank >= 1 && nx_ > 1) ? nx_ + 2*gx_ + cx_ : 1;
  if (my) (*my) = (rank >= 2 && ny_ > 1) ? ny_ + 2*gy_ + cy_ : 1;
  if (mz) (*mz) = (rank >= 3 && nz_ > 1) ? nz_ + 2*gz_ + cz_ : 1;
}
   
//----------------------------------------------------------------------

void FaceFluxes::set_fluxes ( std::vector<double> array, int dx, int dy, int dz)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  ASSERT5("FaceFluxes::set_fluxes",
          "Input array size %lu is smaller than required size %d = %d * %d *%d",
          array.size(),mx*my*mz,mx,my,mz,
          (array.size() >= mx*my*mz) );
  ASSERT4("FaceFluxes::set_fluxes",
          "Input array size %lu is too small for axis strides (dx,dy,dz)=(%d,%d,%d)",
          array.size(),dx,dy,dz,
          (array.size() >= (mx-1)*dx+(my-1)*dy+(mz-1)*dz) );

  ASSERT5("FaceFluxes::set_fluxes",
          "Flux array size %lu is smaller than required size %d = %d * %d *%d",
          fluxes_.size(),mx*my*mz,mx,my,mz,
          (fluxes_.size() >= mx*my*mz) );

    
  for (int ix=0; ix<mx; ix++) {
    for (int iy=0; iy<my; iy++) {
      for (int iz=0; iz<mz; iz++) {
        int i_f = ix + mx*(iy+my*iz);
        int i_a = ix*dx + iy*dy + iz*dz;
        fluxes_[i_f] = array[i_a];
      }
    }
  }
}
  
//----------------------------------------------------------------------

std::vector<double> & FaceFluxes::get_fluxes (int * dx, int * dy, int *dz)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  if (dx) (*dx) = 1;
  if (dy) (*dy) = mx;
  if (dz) (*dz) = mx*my;

  return fluxes_;
}
  
//----------------------------------------------------------------------

float ratio_cell_width (const FaceFluxes & ff_1, const FaceFluxes & ff_2)
{
  double hx_1,hy_1,hz_1;
  double hx_2,hy_2,hz_2;

  ff_1.get_element_size(&hx_1,&hy_1,&hz_1);
  ff_2.get_element_size(&hx_2,&hy_2,&hz_2);
  
  const int rank = ff_1.face().rank();
  ASSERT6("ratio_time_step (FaceFluxes,FaceFluxes)",
          "Anisotropic refinement is not supported %g %g %g  %g %g %g",
          hx_1,hy_1,hz_1,hx_2,hy_2,hz_2,
          ((rank < 2) ||
           ((rank < 3) && (hx_1/hx_2 == hy_1/hy_2)) ||
           ((rank >= 3) && (hx_1/hx_2 == hz_1/hz_2))));
  return hx_1 / hx_2;
}

//----------------------------------------------------------------------

float ratio_time_step (const FaceFluxes & ff_1, const FaceFluxes & ff_2)
{
  return ff_1.time_step() / ff_2.time_step();
}
  
//----------------------------------------------------------------------

void FaceFluxes::coarsen ()
{
  std::vector<double> fluxes_fine = fluxes_;
  int mxf,myf,mzf;
  get_dimensions(&mxf,&myf,&mzf);
  // not correct??  need nxf/2 for ghost and centering
  const int rank = face_.rank();
  // nxc + a = nx/2 + a
  
  const int mxc = (rank >= 1 && nx_ > 1) ? (nx_/2 + 2*gx_ + cx_) : 1;
  const int myc = (rank >= 1 && ny_ > 1) ? (ny_/2 + 2*gy_ + cy_) : 1;
  const int mzc = (rank >= 1 && nz_ > 1) ? (nz_/2 + 2*gz_ + cz_) : 1;

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
#ifdef DEBUG_FACE_FLUXES        
        CkPrintf ("coarsen x-face fluxes_c[%d %d %d %d] += fluxes_f[%d] \n",
                  iycm*dyc + izcm*dzc,
                  iycm*dyc + izcp*dzc,
                  iycp*dyc + izcm*dzc,
                  iycp*dyc + izcp*dzc,i_f);
#endif        
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
#ifdef DEBUG_FACE_FLUXES        
        CkPrintf ("coarsen y-face fluxes_c[%d %d %d %d] += fluxes_f[%d] \n",
                  ixcm*dxc + izcm*dzc,
                  ixcm*dxc + izcp*dzc,
                  ixcp*dxc + izcm*dzc,
                  ixcp*dxc + izcp*dzc,i_f);
#endif        
      }
    }
  } else if (mzf == 1) {
    int izf=0;
    for (int ixf=0; ixf<mxf; ixf++) {
      int ixcm = ixf>>1;
      int ixcp = (ixf+cx_)>>1;
      for (int iyf=0; iyf<myf; iyf++) {
        int iycm = iyf>>1;
        int iycp = (iyf+cy_)>>1;
        int i_f = ixf + mxf*(iyf + myf*izf);
        fluxes_[ixcm*dxc + iycm*dyc] += fluxes_fine[i_f];
        fluxes_[ixcp*dxc + iycm*dyc] += fluxes_fine[i_f];
        fluxes_[ixcm*dxc + iycp*dyc] += fluxes_fine[i_f];
        fluxes_[ixcp*dxc + iycp*dyc] += fluxes_fine[i_f];
#ifdef DEBUG_FACE_FLUXES        
        CkPrintf ("coarsen z-face fluxes_c[%d %d %d %d] += fluxes_f[%d] \n",
                  iycm*dyc + ixcm*dxc,
                  iycm*dyc + ixcp*dxc,
                  iycp*dyc + ixcm*dxc,
                  iycp*dyc + ixcp*dxc,i_f);
#endif
      }
    }
  }
  nx_ = (nx_ != 1) ? nx_>>1 : 1;
  ny_ = (ny_ != 1) ? ny_>>1 : 1;
  nz_ = (nz_ != 1) ? nz_>>1 : 1;
  hx_ = 2.0*hx_;
  hy_ = 2.0*hy_;
  hz_ = 2.0*hz_;
}
  
//----------------------------------------------------------------------

FaceFluxes & FaceFluxes::operator += (const FaceFluxes & face_fluxes_1)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  int mx1,my1,mz1;
  face_fluxes_1.get_dimensions(&mx1,&my1,&mz1);

  ASSERT ("FaceFluxes::operator += (FaceFluxes)",
          "Faces must agree",
          (face() == face_fluxes_1.face()));
  // ASSERT ("FaceFluxes::operator += (FaceFluxes)",
  //         "Sizes must agree",
  //         (mx==mx1) && (my==my1) && (mz==mz1));

  int ixl,ixu,iyl,iyu,izl,izu;
  get_limits(&ixl,&ixu,&iyl,&iyu,&izl,&izu);

  for (int iz=izl; iz<izu; iz++) {
    for (int iy=iyl; iy<iyu; iy++) {
      for (int ix=ixl; ix<ixu; ix++) {
        int i=ix+mx*(iy+my*iz);
        fluxes_[i] += face_fluxes_1.fluxes_[i];
      }
    }
  }
  
  return *this;
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
  
//----------------------------------------------------------------------

FaceFluxes operator - (const FaceFluxes & ff_1, const FaceFluxes & ff_2)
{
  INCOMPLETE("FaceFluxes::operator - (FaceFluxes,FaceFluxes)");
  return ff_1;
}

