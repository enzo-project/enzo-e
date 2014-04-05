// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ScalarExpr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the ScalarExpr class

#ifndef PROBLEM_SCALAR_EXPR_HPP
#define PROBLEM_SCALAR_EXPR_HPP

class ScalarExpr {

  /// @class    ScalarExpr
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  ScalarExpr() throw()
  : param_(0),
    value_(0)
  { }

  /// Destructor
  ~ScalarExpr() throw()
  { }

  /// Copy constructor
  ScalarExpr(const ScalarExpr & scalar_expr) throw()
  { copy_(scalar_expr); }

  /// Assignment operator
  ScalarExpr & operator= (const ScalarExpr & scalar_expr) throw()
  { copy_(scalar_expr); return *this; }

  /// Clone the object
  virtual ScalarExpr * clone() const
  { return (new ScalarExpr(*this)); }
  
  ScalarExpr(Param * param) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    WARNING("MaskExpr::pup()","UNFINISHED");
    p | value_;
    // NOTE: change this function whenever attributes change
  }

  /// Evaluate mask at a point
  double evaluate (double t, double x, double y, double z,
		   Mask * mask = 0, double deflt = 0) const;

  /// Return mask values in an array
  template <class T>
  void evaluate (T *value, double t,
		 int ndx, int nx, double * x,
		 int ndy, int ny, double * y,
		 int ndz, int nz, double * z, 
		 Mask * mask = 0, T * deflt = 0) const;

  /// Return mask values in an array
  template <class T>
  void evaluate (T *value, double t,
		 int ndx, int nx, double * x,
		 int ndy, int ny, double * y,
		 int ndz, int nz, double * z) const
  {
    evaluate(value,t,ndx,nx,x,ndy,ny,y,ndz,nz,z,0,0);
  }

  
private: // functions

  void copy_(const ScalarExpr & scalar_expr) throw();

private: // attributes

  /// Value if param_ type is precision_float_expr
  Param * param_;

  /// Value if param_ type is precision_float
  double value_;

};

template <class T>
void ScalarExpr::evaluate (T * value, double t,
			   int ndx, int nx, double * xv,
			   int ndy, int ny, double * yv,
			   int ndz, int nz, double * zv, 
			   Mask * mask, T * deflt) const
{
  ASSERT6("ScalarExpr::evaluate",
	  "value dimension (%d %d %d) needs to be at least (%d %d %d)",
	  ndx,ndy,ndz,nx,ny,nz,
	  (ndx >= nx) && (ndy >= ny) && (ndz >= nz));

  double * x = new double [ nx*ny*nz ];
  double * y = new double [ nx*ny*nz ];
  double * z = new double [ nx*ny*nz ];
  double * value_temp = new double [nx*ny*nz];
  bool * mv = 0;
  if (mask) {
    mv = new bool [ nx*ny*nz ];
    mask->evaluate(mv, t, nx,nx,xv, ny,ny,yv, nz,nz,zv);
  }
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	x[i] = xv[ix];
	y[i] = yv[iy];
	z[i] = zv[iz];
	value_temp[i] = 0.0;
      }
    }
  }

  int n=nx*ny*nz;

  if (param_) {
    param_->evaluate_float(n, value_temp, x,y,z,t);
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

  if (mv) delete [] mv; mv = 0;
  delete [] value_temp;
  delete [] z;
  delete [] y;
  delete [] x;

}

#endif /* PROBLEM_SCALAR_EXPR_HPP */

