// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixLaplace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     2015-04-02
/// @brief    Implementation of the discrete Laplace operator EnzoMatrixLaplace

#include "enzo.hpp"

// #define DEBUG_MATRIX

//======================================================================

void EnzoMatrixLaplace::matvec (int id_y, int id_x, Block * block,
				int g0) throw()
{
  Data * data = block->data();
  Field field = data->field();

  field.dimensions(0,&mx_,&my_,&mz_);
  data->field_cell_width(&hx_,&hy_,&hz_);

  rank_ = block->rank();

  int precision = field.precision(0);

  void * X = field.values(id_x);
  void * Y = field.values(id_y);

  if      (precision == precision_single)    
    matvec_((float *)(Y),(float *)(X),g0);
  else if (precision == precision_double)    
    matvec_((double *)(Y),(double *)(X),g0);
  else if (precision == precision_quadruple) 
    matvec_((long double *)(Y),(long double *)(X),g0);
  else 
    ERROR1("EnzoMatrixLaplace::matvec()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

void EnzoMatrixLaplace::diagonal (int id_x, Block * block, int g0) throw()
{
  Data * data = block->data();
  Field field = data->field();

  field.dimensions (id_x,&mx_,&my_,&mz_);
  data->field_cell_width(&hx_,&hy_,&hz_);

  rank_ = block->rank();

  int precision = field.precision(0);

  void * X = field.values(id_x);

  if      (precision == precision_single)    
    diagonal_((float *)(X),g0);
  else if (precision == precision_double)    
    diagonal_((double *)(X),g0);
  else if (precision == precision_quadruple) 
    diagonal_((long double *)(X),g0);
  else 
    ERROR1("EnzoMatrixLaplace::diagonal()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMatrixLaplace::matvec_ (T * Y, T * X, int g0) const throw()
{
  const int idx = 1;
  const int idy = mx_;
  const int idz = mx_*my_;
  const int idx2 = 2*idx;
  const int idy2 = 2*idy;
  const int idz2 = 2*idz;

  if (order_ == 2) {

    g0 = std::max(1,g0);

    if (rank_ == 1) {

      for (int ix=g0; ix<mx_-g0; ix++) {
	const int i = ix;
	Y[i] = ( X[i-idx] - 2.0*X[i] + X[i+idx] ) / (hx_*hx_);
      }

    } else if (rank_ == 2) {
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  const int i = ix + mx_*iy;
	  Y[i] = ( X[i+idx] - 2.0*X[i] + X[i-idx]) / (hx_*hx_)
	    +    ( X[i+idy] - 2.0*X[i] + X[i-idy]) / (hy_*hy_);
	}
      }

    } else if (rank_ == 3) {

      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    const int i = ix + mx_*(iy + my_*iz);
	    Y[i] = ( X[i+idx] - 2.0*X[i] + X[i-idx]) / (hx_*hx_)
	      +    ( X[i+idy] - 2.0*X[i] + X[i-idy]) / (hy_*hy_)
	      +    ( X[i+idz] - 2.0*X[i] + X[i-idz]) / (hz_*hz_);
	  }
	}
      }
    }

  } else if (order_ == 4) {

    g0 = std::max(2,g0);

    if (rank_ == 1) {

      const T c0 = -30.0;
      const T c1 = 16.0;
      const T c2 = -1.0;
      const T d  = 12.0*hx_*hx_;
      for (int ix=g0; ix<mx_-g0; ix++) {
	const int i = ix;
	Y[i] = (c0*(X[i]) +
		c1*(X[i-idx] +X[i+idx]) +
		c2*(X[i-idx2]+X[i+idx2])) / d;
      }

    } else if (rank_ == 2) {
      const T c0 = -30.0;
      const T c1 = 16.0;
      const T c2 = -1.0;
      const T dx  = 12.0*hx_*hx_;
      const T dy  = 12.0*hy_*hy_;
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  const int i = ix + mx_*iy;
	  Y[i] = (c0*(X[i]) +
		  c1*(X[i-idx] +X[i+idx]) +
		  c2*(X[i-idx2]+X[i+idx2])) / dx
	    +    (c0*(X[i]) +
		  c1*(X[i-idy] +X[i+idy]) +
		  c2*(X[i-idy2]+X[i+idy2])) / dy;
	}
      }

    } else if (rank_ == 3) {

      const T c0 = -30.0;
      const T c1 = 16.0;
      const T c2 = -1.0;
      const T dx  = 12.0*hx_*hx_;
      const T dy  = 12.0*hy_*hy_;
      const T dz  = 12.0*hz_*hz_;
      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    const int i = ix + mx_*(iy + my_*iz);
	    Y[i] = (c0*(X[i]) +
		    c1*(X[i-idx] +X[i+idx]) +
		    c2*(X[i-idx2]+X[i+idx2])) / dx
	      +    (c0*(X[i]) +
		    c1*(X[i-idy] +X[i+idy]) +
		    c2*(X[i-idy2]+X[i+idy2])) / dy
	      +    (c0*(X[i]) +
		    c1*(X[i-idz] +X[i+idz]) +
		    c2*(X[i-idz2]+X[i+idz2])) / dz;
	  }
	}
      }
    }
  } else {
    ERROR1 ("EnzoMatrixLaplace::diagonal()",
	    "Order %d operator is not supported",
	    order_);
  } 
}

//----------------------------------------------------------------------

template <class T>
void EnzoMatrixLaplace::diagonal_ (T * X, int g0) const throw()
{
  if (order_ == 2) {
    
    g0 = std::max(1,g0);

    // Second-order 7-point discretization
    
    if (rank_ == 1) {
      for (int ix=g0; ix<mx_-g0; ix++) {
	int i = ix;
	X[i] = - 2.0 / (hx_*hx_);
      }
    } else if (rank_ == 2) {
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  int i = ix + mx_*iy;
	  X[i] = - 2.0 / (hx_*hx_)
	    - 2.0 / (hy_*hy_);
	}
      }
    } else if (rank_ == 3) {
      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    X[i] = - 2.0 / (hx_*hx_)
	      +    - 2.0 / (hy_*hy_)
	      +    - 2.0 / (hz_*hz_);
	  }
	}
      }
    }

  } else if (order_ == 4) {

    g0 = std::max(2,g0);

    // Fourth-order 13-point discretization

    if (rank_ == 1) {
      const T c0 = -30.0;
      const T dx  = 12.0*hx_*hx_;
      for (int ix=g0; ix<mx_-g0; ix++) {
	int i = ix;
	X[i] = c0 / dx;
      }
    } else if (rank_ == 2) {
      const T c0 = -30.0;
      const T dx  = 12.0*hx_*hx_;
      const T dy  = 12.0*hy_*hy_;
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  int i = ix + mx_*iy;
	  X[i] = c0 / dx
	    +    c0 / dy;
	}
      }
    } else if (rank_ == 3) {
      const T c0 = -30.0;
      const T dx  = 12.0*hx_*hx_;
      const T dy  = 12.0*hy_*hy_;
      const T dz  = 12.0*hz_*hz_;
      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    X[i] = c0 / dx
	      +    c0 / dy
	      +    c0 / dz;
	  }
	}
      }
    }
  } else {
    ERROR1 ("EnzoMatrixLaplace::diagonal()",
	    "Order %d operator is not supported",
	    order_);
  }
  
}

