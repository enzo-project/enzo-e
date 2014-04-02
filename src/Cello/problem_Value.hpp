// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Value.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the Value class

#ifndef PROBLEM_VALUE_HPP
#define PROBLEM_VALUE_HPP

class Value {

  /// @class    Value
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  Value() throw()
  : mask_list_() 
  { };

  /// Destructor
  ~Value() throw() { };

  /// Copy constructor
  Value(const Value & value) throw()
  { copy_(value); }

  /// Assignment operator
  Value & operator= (const Value & value) throw()
  { copy_(value); return *this; }

  Value(Parameters * parameters,
	const std::string parameter_name) throw();


  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }

  void evaluate
  (double * values, 
   double t,
   int ndx, int nx, double * x,
   int ndy, int ny, double * y,
   int ndz, int nz, double * z) throw ();

  double evaluate (double t, double x, double y, double z) throw ();

private: // functions

  void copy_(const Value & value) throw();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// List of scalar expressions
  std::vector<ScalarExpr *> scalar_expr_list_;

  /// List of masks for each scalar expression
  std::vector<Mask *> mask_list_;

};

#endif /* PROBLEM_VALUE_HPP */

