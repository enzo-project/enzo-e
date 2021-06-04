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
  : scalar_expr_list_(), mask_list_()
  { };

  Value(Parameters * parameters,
	const std::string parameter_name) throw();

  /// Move constructor
  Value( Value &&other) throw() : Value() {swap(*this,other);}

  /// Move assignment. Replace the contents of *this with that of other
  Value & operator=(Value &&other) throw() { swap(*this,other); return *this;}

  /// delete the copy constructor and copy assignment operator
  Value(const Value & value) = delete;
  Value & operator= (const Value & value) = delete;

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  template <class T>
  void evaluate
  (T * values, double t,
   int ndx, int nx, double * x,
   int ndy, int ny, double * y,
   int ndz, int nz, double * z) throw ();

  double evaluate (double t, double x, double y, double z) throw ();

  /// Swaps the contents of the first Value object with another
  friend void swap(Value &first, Value &second){
    // NOTE: change this function whenever attributes change
    std::swap(first.scalar_expr_list_, second.scalar_expr_list_);
    std::swap(first.mask_list_, second.mask_list_);
  }

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// List of scalar expressions
  std::vector<ScalarExpr> scalar_expr_list_;

  /// List of masks for each scalar expression
  ///
  /// Note: Entries are only ever shared temporarily.
  std::vector< std::shared_ptr<Mask> > mask_list_;

};

#endif /* PROBLEM_VALUE_HPP */
