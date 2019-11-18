// See LICENSE_CELLO file for license and copyright information

/// @file     data_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    [\ref Data] Implementation of the FaceFluxes class

#include "data.hpp"

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
  face_.get_dimensions(&nx_,&ny_,&nz_,nx,ny,nz);
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
  if (nx_ == 1) gx_ = 0;
  if (ny_ == 1) gy_ = 0;
  if (nz_ == 1) gz_ = 0;
}

//----------------------------------------------------------------------

void FaceFluxes::set_centering(int cx, int cy, int cz)
{
  cx_ = cx;
  cy_ = cy;
  cz_ = cz;
  if (nx_ == 1) cx_ = 0;
  if (ny_ == 1) cy_ = 0;
  if (nz_ == 1) cz_ = 0;
}
  
//----------------------------------------------------------------------

void FaceFluxes::allocate ()
{
  int mx,my,mz;
  get_dimensions (&mx,&my,&mz);
  if (mx*my*mz>0) fluxes_.resize(mx*my*mz);
}

//----------------------------------------------------------------------

void FaceFluxes::deallocate()
{ }

//----------------------------------------------------------------------

Face FaceFluxes::face () const
{ return face_; }

//----------------------------------------------------------------------

void FaceFluxes::get_element_size (double *hx, double *hy, double * hz) const
{
  if (hx) (*hx) = hx_;
  if (hy) (*hy) = hy_;
  if (hz) (*hz) = hz_;
}
  
//----------------------------------------------------------------------

double FaceFluxes::time_step () const
{ return dt_;}
  
//----------------------------------------------------------------------

void FaceFluxes::get_dimensions (int *mx, int *my, int *mz) const
{
  if (mx) (*mx) = nx_ + 2*gx_ + cx_;
  if (my) (*my) = ny_ + 2*gy_ + cy_;
  if (mz) (*mz) = nz_ + 2*gz_ + cz_;
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
{ return 0.0;}

//----------------------------------------------------------------------

float ratio_time_step (const FaceFluxes & ff_1, const FaceFluxes & ff_2)
{ return 0.0; }
  
//----------------------------------------------------------------------

void FaceFluxes::coarsen ()
{ }
  
//----------------------------------------------------------------------

FaceFluxes & FaceFluxes::operator += (const FaceFluxes & face_fluxes)
{ return *this; }
  
//----------------------------------------------------------------------

FaceFluxes & FaceFluxes::operator *= (double weight)
{ return *this; }
  
//----------------------------------------------------------------------

FaceFluxes operator - (const FaceFluxes & ff_1, const FaceFluxes & ff_2)
{ return ff_1; }

