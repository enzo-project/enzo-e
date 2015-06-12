// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixLaplace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-04-02
/// @brief    Implementation of the discrete Laplace operator EnzoMatrixLaplace

#include "enzo.hpp"

//======================================================================

void EnzoMatrixLaplace::matvec (int id_y, int id_x, Block * block) throw()
{
  if (! block->is_leaf()) return;

  Data * data = block->data();
  Field field = data->field();

  field.dimensions(0,&mx_,&my_,&mz_);
  field.size        (&nx_,&ny_,&nz_);
  field.ghost_depth    (0,&gx_,&gy_,&gz_);
  data->field_cell_width(&hx_,&hy_,&hz_);

  rank_ = (my_ == 1) ? 1 : (mz_ == 1) ? 2 : 3;

  int precision = field.precision(0);

  void * X = field.values(id_x);
  void * Y = field.values(id_y);

  if      (precision == precision_single)    
    matvec_((float *)(Y),(float *)(X));
  else if (precision == precision_double)    
    matvec_((double *)(Y),(double *)(X));
  else if (precision == precision_quadruple) 
    matvec_((long double *)(Y),(long double *)(X));
  else 
    ERROR1("EnzoMethodGravityCg()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

void EnzoMatrixLaplace::diagonal (int id_x, Block * block) throw()
{
  if (! block->is_leaf()) return;

  Data * data = block->data();
  Field field = data->field();

  field.dimensions (id_x,&mx_,&my_,&mz_);
  field.size            (&nx_,&ny_,&nz_);
  field.ghost_depth     (id_x,&gx_,&gy_,&gz_);
  data->field_cell_width(&hx_,&hy_,&hz_);

  rank_ = (my_ == 1) ? 1 : (mz_ == 1) ? 2 : 3;

  int precision = field.precision(0);

  void * X = field.values(id_x);

  if      (precision == precision_single)    
    diagonal_((float *)(X));
  else if (precision == precision_double)    
    diagonal_((double *)(X));
  else if (precision == precision_quadruple) 
    diagonal_((long double *)(X));
  else 
    ERROR1("EnzoMethodGravityCg()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMatrixLaplace::matvec_ (T * Y, T * X) const throw()
{
  const int idx = 1;
  const int idy = mx_;
  const int idz = mx_*my_;

  const int i0 = gx_ + mx_*(gy_ + my_*gz_);

  if (rank_ == 1) {
    for (int ix=0; ix<nx_; ix++) {
      int i = i0 + ix;
      Y[i] = ( X[i-idx] - 2.0*X[i] + X[i+idx] ) / (hx_*hx_);
    }
  } else if (rank_ == 2) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + ix + mx_*iy;
	Y[i] = ( X[i+idx] - 2.0*X[i] + X[i-idx]) / (hx_*hx_)
	  +    ( X[i+idy] - 2.0*X[i] + X[i-idy]) / (hy_*hy_);
      }
    }
  } else if (rank_ == 3) {
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + ix + mx_*(iy + my_*iz);
	  Y[i] = ( X[i+idx] - 2.0*X[i] + X[i-idx]) / (hx_*hx_)
	    +    ( X[i+idy] - 2.0*X[i] + X[i-idy]) / (hy_*hy_)
	    +    ( X[i+idz] - 2.0*X[i] + X[i-idz]) / (hz_*hz_);
	}
      }
    }
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoMatrixLaplace::diagonal_ (T * X) const throw()
{
  const int i0 = gx_ + mx_*(gy_ + my_*gz_);

  if (rank_ == 1) {
    for (int ix=0; ix<nx_; ix++) {
      int i = i0 + ix;
      X[i] = - 2.0 / (hx_*hx_);
    }
  } else if (rank_ == 2) {
    for (int iy=0; iy<ny_; iy++) {
      for (int ix=0; ix<nx_; ix++) {
	int i = i0 + ix + mx_*iy;
	X[i] = - 2.0 / (hx_*hx_)
	       - 2.0 / (hy_*hy_);
      }
    }
  } else if (rank_ == 3) {
    for (int iz=0; iz<nz_; iz++) {
      for (int iy=0; iy<ny_; iy++) {
	for (int ix=0; ix<nx_; ix++) {
	  int i = i0 + ix + mx_*(iy + my_*iz);
	  X[i] = - 2.0 / (hx_*hx_)
	    +    - 2.0 / (hy_*hy_)
	    +    - 2.0 / (hz_*hz_);
	}
      }
    }
  }
}

