// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixLaplace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     2015-04-02
/// @brief    Implementation of the discrete Laplace operator EnzoMatrixLaplace

#include "enzo.hpp"

// #define DEBUG_MATRIX

//======================================================================

void EnzoMatrixLaplace::matvec (int i_y, int i_x, Block * block,
				int g0) throw()
{
  Field field = block->data()->field();

  field.dimensions(0,&mx_,&my_,&mz_);
  block->cell_width (&hx_,&hy_,&hz_);
  
  rank_ = block->rank();

  enzo_float * X = (enzo_float * ) field.values(i_x);
  enzo_float * Y = (enzo_float * ) field.values(i_y);
  
  matvec_(Y,X,g0);
}

//----------------------------------------------------------------------

void EnzoMatrixLaplace::diagonal (int i_x, Block * block, int g0) throw()
{
  Field field = block->data()->field();

  field.dimensions (i_x,&mx_,&my_,&mz_);
  block->cell_width     (&hx_,&hy_,&hz_);

  rank_ = block->rank();

  enzo_float * X = (enzo_float * ) field.values(i_x);

  diagonal_(X,g0);
}

//----------------------------------------------------------------------

void EnzoMatrixLaplace::matvec_ (enzo_float * Y, enzo_float * X, int g0) const throw()
{
  const int idx = 1;
  const int idy = mx_;
  const int idz = mx_*my_;

  if (order_ == 2) {

    g0 = std::max(1,g0);

    double dx = (rank_ >= 1) ? 1.0 / (hx_*hx_) : 0.0;
    double dy = (rank_ >= 2) ? 1.0 / (hy_*hy_) : 0.0;
    double dz = (rank_ >= 3) ? 1.0 / (hz_*hz_) : 0.0;

    if (rank_ == 1) {
      for (int ix=g0; ix<mx_-g0; ix++) {
	const int i = ix;
	Y[i] = ( X[i-idx] - 2.0*X[i] + X[i+idx] ) * dx;
      }

    } else if (rank_ == 2) {
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  const int i = ix + mx_*iy;
	  Y[i] = ( X[i+idx] - 2.0*X[i] + X[i-idx]) * dx
	    +    ( X[i+idy] - 2.0*X[i] + X[i-idy]) * dy;
	}
      }

    } else if (rank_ == 3) {
      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    const int i = ix + mx_*(iy + my_*iz);
	    Y[i] = ( X[i+idx] - 2.0*X[i] + X[i-idx]) * dx
	      +    ( X[i+idy] - 2.0*X[i] + X[i-idy]) * dy
	      +    ( X[i+idz] - 2.0*X[i] + X[i-idz]) * dz;
	  }
	}
      }
    }

  } else if (order_ == 4) {

    const int idx2 = 2*idx;
    const int idy2 = 2*idy;
    const int idz2 = 2*idz;

    g0 = std::max(2,g0);

    const enzo_float c0 = -30.0;
    const enzo_float c1 = 16.0;
    const enzo_float c2 = -1.0;
    const enzo_float dx = (rank_ >= 1) ? 1.0/(12.0*hx_*hx_) : 0.0;
    const enzo_float dy = (rank_ >= 2) ? 1.0/(12.0*hy_*hy_) : 0.0;
    const enzo_float dz = (rank_ >= 3) ? 1.0/(12.0*hz_*hz_) : 0.0;

    if (rank_ == 1) {

      for (int ix=g0; ix<mx_-g0; ix++) {
	const int i = ix;
	Y[i] = (c0*(X[i]) +
		c1*(X[i-idx] +X[i+idx]) +
		c2*(X[i-idx2]+X[i+idx2])) * dx;
      }

    } else if (rank_ == 2) {

      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  const int i = ix + mx_*iy;
	  Y[i] = (c0*(X[i]) +
		  c1*(X[i-idx] +X[i+idx]) +
		  c2*(X[i-idx2]+X[i+idx2])) * dx
	    +    (c0*(X[i]) +
		  c1*(X[i-idy] +X[i+idy]) +
		  c2*(X[i-idy2]+X[i+idy2])) * dy;
	}
      }

    } else if (rank_ == 3) {

      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    const int i = ix + mx_*(iy + my_*iz);
	    Y[i] = (c0*(X[i]) +
		    c1*(X[i-idx] +X[i+idx]) +
		    c2*(X[i-idx2]+X[i+idx2])) * dx
	      +    (c0*(X[i]) +
		    c1*(X[i-idy] +X[i+idy]) +
		    c2*(X[i-idy2]+X[i+idy2])) * dy
	      +    (c0*(X[i]) +
		    c1*(X[i-idz] +X[i+idz]) +
		    c2*(X[i-idz2]+X[i+idz2])) * dz;
	  }
	}
      }
    }
  } else if (order_ == 6) {

    const int idx2 = 2*idx;
    const int idy2 = 2*idy;
    const int idz2 = 2*idz;
    const int idx3 = 3*idx;
    const int idy3 = 3*idy;
    const int idz3 = 3*idz;

    g0 = std::max(3,g0);

    const enzo_float c0 = -2720.0;
    const enzo_float c1 = 1455.0;
    const enzo_float c2 = -96.0;
    const enzo_float c3 = 1.0;
    const enzo_float dx = (rank_ >= 1) ? 1.0/(1080.0*hx_*hx_) : 0.0;
    const enzo_float dy = (rank_ >= 2) ? 1.0/(1080.0*hy_*hy_) : 0.0;
    const enzo_float dz = (rank_ >= 3) ? 1.0/(1080.0*hz_*hz_) : 0.0;

    if (rank_ == 1) {

      for (int ix=g0; ix<mx_-g0; ix++) {
	const int i = ix;
	Y[i] = (c0*(X[i]) +
		c1*(X[i-idx] +X[i+idx]) +
		c2*(X[i-idx2]+X[i+idx2]) +
		c3*(X[i-idx3]+X[i+idx3])) * dx;
      }

    } else if (rank_ == 2) {

      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  const int i = ix + mx_*iy;
	  Y[i] = (c0*(X[i]) +
		  c1*(X[i-idx] +X[i+idx]) +
		  c2*(X[i-idx2]+X[i+idx2]) +
		  c3*(X[i-idx3]+X[i+idx3])) * dx
	    +    (c0*(X[i]) +
		  c1*(X[i-idy] +X[i+idy]) +
		  c2*(X[i-idy2]+X[i+idy2]) +
		  c3*(X[i-idy3]+X[i+idy3])) * dy;
	}
      }

    } else if (rank_ == 3) {

      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    const int i = ix + mx_*(iy + my_*iz);
	    Y[i] = (c0*(X[i]) +
		    c1*(X[i-idx] +X[i+idx]) +
		    c2*(X[i-idx2]+X[i+idx2]) +
		    c3*(X[i-idx3]+X[i+idx3])) * dx
	      +    (c0*(X[i]) +
		    c1*(X[i-idy] +X[i+idy]) +
		    c2*(X[i-idy2]+X[i+idy2]) +
		    c2*(X[i-idy3]+X[i+idy3])) * dy
	      +    (c0*(X[i]) +
		    c1*(X[i-idz] +X[i+idz]) +
		    c2*(X[i-idz2]+X[i+idz2]) +
		    c3*(X[i-idz3]+X[i+idz3])) * dz;
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

void EnzoMatrixLaplace::diagonal_ (enzo_float * X, int g0) const throw()
{
  if (order_ == 2) {
    
    g0 = std::max(1,g0);

    // Second-order 7-point discretization
    
    double dx = (rank_ >= 1) ? 1.0/(hx_*hx_) : 0.0;
    double dy = (rank_ >= 2) ? 1.0/(hy_*hy_) : 0.0;
    double dz = (rank_ >= 3) ? 1.0/(hz_*hz_) : 0.0;
    
    if (rank_ == 1) {
      for (int ix=g0; ix<mx_-g0; ix++) {
	int i = ix;
	X[i] = - 2.0 * dx;
      }
    } else if (rank_ == 2) {
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  int i = ix + mx_*iy;
	  X[i] = - 2.0 * dx
	    -      2.0 * dy;
	}
      }
    } else if (rank_ == 3) {
      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    X[i] = - 2.0 * dx
	      +    - 2.0 * dy
	      +    - 2.0 * dz;
	  }
	}
      }
    }

  } else if (order_ == 4) {

    const enzo_float c0 = -30.0;
    const enzo_float dx = (rank_ >= 1) ? 1.0/(12.0*hx_*hx_) : 0.0;
    const enzo_float dy = (rank_ >= 2) ? 1.0/(12.0*hy_*hy_) : 0.0;
    const enzo_float dz = (rank_ >= 3) ? 1.0/(12.0*hz_*hz_) : 0.0;

    g0 = std::max(2,g0);

    // Fourth-order 13-point discretization

    if (rank_ == 1) {
      for (int ix=g0; ix<mx_-g0; ix++) {
	int i = ix;
	X[i] = c0 * dx;
      }
    } else if (rank_ == 2) {
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  int i = ix + mx_*iy;
	  X[i] = c0 * dx
	    +    c0 * dy;
	}
      }
    } else if (rank_ == 3) {
      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    X[i] = c0 * dx
	      +    c0 * dy
	      +    c0 * dz;
	  }
	}
      }
    }

  } else if (order_ == 6) {

    g0 = std::max(3,g0);

    // Sixth-order 19-point discretization

    const enzo_float c0 = -2720.0;
    const enzo_float dx = (rank_ >= 1) ? 1.0/(1080.0*hx_*hx_) : 0.0;
    const enzo_float dy = (rank_ >= 2) ? 1.0/(1080.0*hy_*hy_) : 0.0;
    const enzo_float dz = (rank_ >= 3) ? 1.0/(1080.0*hz_*hz_) : 0.0;
      
    if (rank_ == 1) {
      
      for (int ix=g0; ix<mx_-g0; ix++) {
	int i = ix;
	X[i] = c0 * dx;
      }
    } else if (rank_ == 2) {
      for   (int iy=g0; iy<my_-g0; iy++) {
	for (int ix=g0; ix<mx_-g0; ix++) {
	  int i = ix + mx_*iy;
	  X[i] = c0 * dx
	    +    c0 * dy;
	}
      }
    } else if (rank_ == 3) {
      for     (int iz=g0; iz<mz_-g0; iz++) {
	for   (int iy=g0; iy<my_-g0; iy++) {
	  for (int ix=g0; ix<mx_-g0; ix++) {
	    int i = ix + mx_*(iy + my_*iz);
	    X[i] = c0 * dx
	      +    c0 * dy
	      +    c0 * dz;
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

