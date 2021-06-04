// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ScalarExpr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the ScalarExpr class

#ifndef PROBLEM_SCALAR_EXPR_HPP
#define PROBLEM_SCALAR_EXPR_HPP

#include "_parameters.hpp" // require full declaration of Expression class

class ScalarExpr {

  /// @class    ScalarExpr
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Default Constructor
  ScalarExpr() throw()
  : expr_(),
    value_(0)
  { }

  /// Destructor
  ~ScalarExpr() = default;

  /// Copy constructor
  ScalarExpr(const ScalarExpr & scalar_expr) = default;

  /// Copy assignment operator
  ScalarExpr & operator= (const ScalarExpr & scalar_expr) = default;

  /// Main Constructor
  ScalarExpr(Parameters * parameters,
	     const std::string &parameter_name,
	     int parameter_index = -1) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | expr_;
    p | value_;
  }

  /// Evaluate mask at a point
  double evaluate (double t, double x, double y, double z,
		   std::shared_ptr<Mask> mask = nullptr, double deflt = 0) const;

  /// Return mask values in an array
  template <class T>
  void evaluate (T *value, double t,
		 int ndx, int nx, double * x,
		 int ndy, int ny, double * y,
		 int ndz, int nz, double * z, 
		 std::shared_ptr<Mask> mask = nullptr, T * deflt = 0) const;

  /// Return mask values in an array
  template <class T>
  void evaluate (T *value, double t,
		 int ndx, int nx, double * x,
		 int ndy, int ny, double * y,
		 int ndz, int nz, double * z) const
  {
    evaluate(value,t,ndx,nx,x,ndy,ny,y,ndz,nz,z,0,0);
  }

private: // attributes

  /// Value if param_ type is precision_float_expr
  Expression expr_;

  /// Value if param_ type is precision_float
  double value_;

};
#endif /* PROBLEM_SCALAR_EXPR_HPP */

