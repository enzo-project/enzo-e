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
  ScalarExpr() throw();

  /// Destructor
  ~ScalarExpr() throw();

  /// Copy constructor
  ScalarExpr(const ScalarExpr & scalar_expr) throw();

  /// Assignment operator
  ScalarExpr & operator= (const ScalarExpr & scalar_expr) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }

  virtual void evaluate
  (FieldBlock * field_block, int index_field) throw ();

  
private: // functions


private: // attributes


};

#endif /* PROBLEM_SCALAR_EXPR_HPP */

