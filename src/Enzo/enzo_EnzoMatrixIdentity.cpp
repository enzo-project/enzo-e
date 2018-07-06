// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixIdentity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     2015-04-02
/// @brief    Implementation of the discrete Identity operator EnzoMatrixIdentity

#include "enzo.hpp"

//======================================================================

void EnzoMatrixIdentity::matvec 
(int id_y, int id_x, Block * block, int g0) throw()
{
  Field field = block->data()->field();

  field.dimensions(0,&mx_,&my_,&mz_);

  enzo_float * X = (enzo_float * ) field.values(id_x);
  enzo_float * Y = (enzo_float * ) field.values(id_y);

  matvec_(Y,X,g0);
}

//----------------------------------------------------------------------

void EnzoMatrixIdentity::diagonal (int id_x, Block * block, int g0) throw()
{
  Field field = block->data()->field();

  field.dimensions(0,&mx_,&my_,&mz_);

  enzo_float * X = (enzo_float * ) field.values(id_x);

  diagonal_ (X,g0);
}

//----------------------------------------------------------------------

void EnzoMatrixIdentity::matvec_ (enzo_float * Y, enzo_float * X, int g0) const throw()
{
  const int ix0 = (mx_ > 1) ? g0 : 0;
  const int iy0 = (my_ > 1) ? g0 : 0;
  const int iz0 = (mz_ > 1) ? g0 : 0;

  for (int iz=iz0; iz<mz_-iz0; iz++) {
    for (int iy=iy0; iy<my_-iy0; iy++) {
      for (int ix=ix0; ix<mx_-ix0; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	Y[i] = X[i];
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMatrixIdentity::diagonal_ (enzo_float * X, int g0) const throw()
{
  const int ix0 = (mx_ > 1) ? g0 : 0;
  const int iy0 = (my_ > 1) ? g0 : 0;
  const int iz0 = (mz_ > 1) ? g0 : 0;

  for (int iz=iz0; iz<mz_-iz0; iz++) {
    for (int iy=iy0; iy<my_-iy0; iy++) {
      for (int ix=ix0; ix<mx_-ix0; ix++) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] = 1.0;
      }
    }
  }
}
