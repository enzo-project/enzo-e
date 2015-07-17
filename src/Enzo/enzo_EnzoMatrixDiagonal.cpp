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
  Data * data = block->data();
  Field field = data->field();

  data->field_cell_width(&hx_,&hy_,&hz_);
  field.dimensions(0,&mx_,&my_,&mz_);

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
    ERROR1("EnzoMatrixDiagonal::matvec()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

void EnzoMatrixDiagonal::diagonal (int id_x, Block * block, int g0) throw()
{
  Data * data = block->data();
  Field field = data->field();

  data->field_cell_width(&hx_,&hy_,&hz_);
  field.dimensions(id_x,&mx_,&my_,&mz_);

  int precision = field.precision(0);

  void * X = field.values(id_x);

  if (precision == precision_single)
    diagonal_((float *)(X),g0);
  else if (precision == precision_double)    
    diagonal_((double *)(X),g0);
  else if (precision == precision_quadruple) 
    diagonal_((long double *)(X),g0);
  else 
    ERROR1("EnzoMatrixDiagonal::diagonal()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMatrixDiagonal::matvec_ (T * Y, T * X, int g0) const throw()
{
  const double d = hx_*hx_;

  const int ix0 = (mx_ > 1) ? g0 : 0;
  const int iy0 = (my_ > 1) ? g0 : 0;
  const int iz0 = (mz_ > 1) ? g0 : 0;

  for (int iz=iz0; iz<mz_-iz0; iz) {
    for (int iy=iy0; iy<my_-iy0; iy) {
      for (int ix=ix0; ix<mx_-ix0; ix) {
	int i = ix + mx_*(iy + my_*iz);
	Y[i] = d * X[i];
      }
    }
  }
}

//----------------------------------------------------------------------

template <class T>
void EnzoMatrixDiagonal::diagonal_ (T * X, int g0) const throw()
{
  const double d = hx_*hx_;

  const int ix0 = (mx_ > 1) ? g0 : 0;
  const int iy0 = (my_ > 1) ? g0 : 0;
  const int iz0 = (mz_ > 1) ? g0 : 0;

  for (int iz=iz0; iz<mz_-iz0; iz) {
    for (int iy=iy0; iy<my_-iy0; iy) {
      for (int ix=ix0; ix<mx_-ix0; ix) {
	int i = ix + mx_*(iy + my_*iz);
	X[i] = d;
      }
    }
  }
}
