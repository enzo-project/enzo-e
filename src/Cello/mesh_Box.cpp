// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-09-21
/// @brief    

#include "mesh.hpp"
//----------------------------------------------------------------------

void Box::get_send_limits (int * ixm, int * ixp,
                           int * iym, int * iyp,
                           int * izm, int * izp)
{

  (*ixm) = (rank_ >= 1) ? ixm_ + gx_ : 0;
  (*ixp) = (rank_ >= 1) ? ixp_ + gx_ : 1;
  (*iym) = (rank_ >= 2) ? iym_ + gy_ : 0;
  (*iyp) = (rank_ >= 2) ? iyp_ + gy_ : 1;
  (*izm) = (rank_ >= 3) ? izm_ + gz_ : 0;
  (*izp) = (rank_ >= 3) ? izp_ + gz_ : 1;
}

//----------------------------------------------------------------------

void Box::get_recv_limits (int * ixm, int * ixp,
                           int * iym, int * iyp,
                           int * izm, int * izp)
{
  (*ixm) = (rank_ >= 1) ? (ixm_ - ix0_)/r_ + gx_ : 0;
  (*ixp) = (rank_ >= 1) ? (ixp_ - ix0_)/r_ + gx_ : 1;
  (*iym) = (rank_ >= 2) ? (iym_ - iy0_)/r_ + gy_ : 0;
  (*iyp) = (rank_ >= 2) ? (iyp_ - iy0_)/r_ + gy_ : 1;
  (*izm) = (rank_ >= 3) ? (izm_ - iz0_)/r_ + gz_ : 0;
  (*izp) = (rank_ >= 3) ? (izp_ - iz0_)/r_ + gz_ : 1;
}

//----------------------------------------------------------------------
  
void Box::compute_region()
                      
{
  if (level_ == -1) {
    ix0_ = nx_*(2.0*fx_ - cx_);
    iy0_ = ny_*(2.0*fy_ - cy_);
    iz0_ = nz_*(2.0*fz_ - cz_);
  } else if (level_ == 0) {
    ix0_ = nx_*fx_;
    iy0_ = ny_*fy_;
    iz0_ = nz_*fz_;
  } else if (level_ == +1) {
    ix0_ = nx_*(fx_ + 0.5*cx_);
    iy0_ = ny_*(fy_ + 0.5*cy_);
    iz0_ = nz_*(fz_ + 0.5*cz_);
  }

  ixm_ = std::max(-gxs_, int(floor(ix0_ - r_*gxr_)));
  iym_ = std::max(-gys_, int(floor(iy0_ - r_*gyr_)));
  izm_ = std::max(-gzs_, int(floor(iz0_ - r_*gzr_)));

  ixp_ = std::min(nx_ + gxs_, int(ceil(ix0_ + r_*(nx_+gxr_))));
  iyp_ = std::min(ny_ + gys_, int(ceil(iy0_ + r_*(ny_+gyr_))));
  izp_ = std::min(nz_ + gzs_, int(ceil(iz0_ + r_*(nz_+gzr_))));
}

