// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ScalarExpr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2014-03-31
/// @brief    Implementation of the ScalarExpr class

#include <cstring>

#include "problem.hpp"

//----------------------------------------------------------------------

ScalarExpr::ScalarExpr (Parameters * parameters,
			const std::string &parameter_name,
			int parameter_index) throw()
  : expr_(),
    value_(0)
{
  Param* param = parameters->param(parameter_name, parameter_index);
  if (param->type() == parameter_float) {
    value_ = param->get_float();
  } else if (parameter_index == -1) {
    expr_ = parameters->value_Expression(parameter_name);
  } else {
    expr_ = parameters->list_value_Expression(parameter_index, parameter_name);
  }
}

//----------------------------------------------------------------------

double ScalarExpr::evaluate (double t, double x, double y, double z, 
			     std::shared_ptr<Mask> mask, double deflt) const
{
  double value;
  bool m = mask ? mask->evaluate(t,x,y,z) : true;
  if (m) {
    if (expr_.initialized())
      expr_.evaluate_float(1,&value,&x,&y,&z,t);
    else
      value = value_;
  } else {
    value = deflt;
  }
  return value;
}

//----------------------------------------------------------------------

template <class T>
void ScalarExpr::evaluate (T * value, double t,
			   int ndx, int nx, double * xv,
			   int ndy, int ny, double * yv,
			   int ndz, int nz, double * zv, 
			   std::shared_ptr<Mask> mask, T * deflt) const
{
  ASSERT6("ScalarExpr::evaluate",
	  "value dimension (%d %d %d) needs to be at least (%d %d %d)",
	  ndx,ndy,ndz,nx,ny,nz,
	  (ndx >= nx) && (ndy >= ny) && (ndz >= nz));

  bool * mv = 0;
  if (mask) {
    mv = new bool [ nx*ny*nz ];
    mask->evaluate(mv, t, nx,nx,xv, ny,ny,yv, nz,nz,zv);
  }

  double * value_temp = new double [nx*ny*nz];
  double * x = new double [ nx*ny*nz ];
  double * y = new double [ nx*ny*nz ];
  double * z = new double [ nx*ny*nz ];

  for (int i=0; i<nx*ny*nz; i++) {
    x[i] = y[i] = z[i] = 0.0;
  }
  
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i=ix + nx*(iy + ny*iz);
	x[i] = xv[ix];
	y[i] = yv[iy];
	z[i] = zv[iz];
	value_temp[i] = 0.0;
      }
    }
  }

  int n=nx*ny*nz;

  if (expr_.initialized()) {
    expr_.evaluate_float(n, value_temp, x,y,z,t);
  } else {
    for (int i=0; i<n; i++) value_temp[i]=value_;
  }

  if (mask) {
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i=ix + nx*(iy + ny*iz);
	  int id=ix + ndx*(iy + ndy*iz);
	  value[id] = mv[i] ? (T) value_temp[i] : deflt[id];
	}
      }
    }
  } else {
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i=ix + nx*(iy + ny*iz);
	  int id=ix + ndx*(iy + ndy*iz);
	  value[id] = (T) (value_temp[i]);
	}
      }
    }
  }

  if (mv) { delete [] mv; }
  delete [] value_temp;
  delete [] z;
  delete [] y;
  delete [] x;

}


template void ScalarExpr::evaluate
(float *value, double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z, 
 std::shared_ptr<Mask> mask, float * deflt) const;
template void ScalarExpr::evaluate
(double *value, double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z, 
 std::shared_ptr<Mask> mask, double * deflt) const;
template void ScalarExpr::evaluate
(long double *value, double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z, 
 std::shared_ptr<Mask> mask, long double * deflt) const;

//----------------------------------------------------------------------

bool ScalarExpr::wraps_single_float_param() const
{ return !expr_.initialized(); }
