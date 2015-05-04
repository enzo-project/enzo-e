// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixDiagonal.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-04-02
/// @brief    Implementation of the discrete Diagonal operator EnzoMatrixDiagonal

#include "enzo.hpp"

//======================================================================

void EnzoMatrixDiagonal::matvec (int id_y, int id_x, Block * block) throw()
{
  if (! block->is_leaf()) return;

  Data * data = block->data();
  Field field = data->field();

  data->field_cell_width(&hx_,&hy_,&hz_);
  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);
  m_ = mx*my*mz;

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

void EnzoMatrixDiagonal::diagonal (int id_x, Block * block) throw()
{
  if (! block->is_leaf()) return;

  Data * data = block->data();
  Field field = data->field();

  data->field_cell_width(&hx_,&hy_,&hz_);
  int mx,my,mz;
  field.dimensions(id_x,&mx,&my,&mz);
  m_ = mx*my*mz;

  int precision = field.precision(0);

  void * X = field.values(id_x);

  if (precision == precision_single)
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
void EnzoMatrixDiagonal::matvec_ (T * Y, T * X) const throw()
{
  const double d = hx_*hx_;

  for (int i=0; i<m_; i++) Y[i] = d * X[i];
}

//----------------------------------------------------------------------

template <class T>
void EnzoMatrixDiagonal::diagonal_ (T * X) const throw()
{
  const double d = hx_*hx_;

  for (int i=0; i<m_; i++) X[i] = d;
}
