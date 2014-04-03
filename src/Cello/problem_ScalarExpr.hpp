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
  void evaluate (double *value, double t,
		 int ndx, int nx, double * x,
		 int ndy, int ny, double * y,
		 int ndz, int nz, double * z, 
		 Mask * mask = 0, double * deflt = 0) const;

  
private: // functions

  void copy_(const ScalarExpr & scalar_expr) throw();

private: // attributes

  /// Value if param_ type is precision_float_expr
  Param * param_;

  /// Value if param_ type is precision_float
  double value_;

};

#endif /* PROBLEM_SCALAR_EXPR_HPP */

