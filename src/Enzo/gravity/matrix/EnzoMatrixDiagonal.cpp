// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixDiagonal.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     2015-04-02
/// @brief    Implementation of the discrete Diagonal operator EnzoMatrixDiagonal

#include "enzo.hpp"

//======================================================================

void EnzoMatrixDiagonal::matvec (int id_y, int id_x, Block * block, int g0) throw()
{
  Field field = block->data()->field();

  block->cell_width (&hx_,&hy_,&hz_);
  field.dimensions(0,&mx_,&my_,&mz_);

  enzo_float * X = (enzo_float * ) field.values(id_x);
  enzo_float * Y = (enzo_float * ) field.values(id_y);

  matvec_(Y,X,g0);
}

//----------------------------------------------------------------------

void EnzoMatrixDiagonal::matvec
(precision_type precision,
 void * y, void * x, int g0) throw()
{
  matvec_((enzo_float *)(y),(enzo_float *)(x),g0);
}

//----------------------------------------------------------------------

void EnzoMatrixDiagonal::diagonal (int id_x, Block * block, int g0) throw()
{
  Field field = block->data()->field();

  block->cell_width    (&hx_,&hy_,&hz_);
  field.dimensions(id_x,&mx_,&my_,&mz_);

  enzo_float * X = (enzo_float * ) field.values(id_x);
  
  diagonal_ (X,g0);
}

//----------------------------------------------------------------------

// template <class T>
void EnzoMatrixDiagonal::matvec_ (enzo_float * Y, enzo_float * X, int g0) const throw()
{
  const double d = hx_*hx_;

  const int ix0 = (mx_ > 1) ? g0 : 0;
  const int iy0 = (my_ > 1) ? g0 : 0;
  const int iz0 = (mz_ > 1) ? g0 : 0;

  for (int iz=iz0; iz<mz_-iz0; iz++) {
    for (int iy=iy0; iy<my_-iy0; iy++) {
      for (int ix=ix0; ix<mx_-ix0; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	Y[i] = d * X[i];
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMatrixDiagonal::diagonal_ (enzo_float * X, int g0) const throw()
{
  const double d = hx_*hx_;

  const int ix0 = (mx_ > 1) ? g0 : 0;
  const int iy0 = (my_ > 1) ? g0 : 0;
  const int iz0 = (mz_ > 1) ? g0 : 0;

  for (int iz=iz0; iz<mz_-iz0; iz++) {
    for (int iy=iy0; iy<my_-iy0; iy++) {
      for (int ix=ix0; ix<mx_-ix0; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] = d;
      }
    }
  }
}
