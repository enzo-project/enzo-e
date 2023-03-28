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

  /// we delete the copy constructor/assignment operator since the ScalarExpr
  /// class doesn't implement it
  Value(const Value & value) = delete;
  Value & operator= (const Value & value) = delete;

  /// Default move constructor and assignment operators
  Value(Value&&) = default;
  Value& operator= (Value&&) = default;

  Value(Parameters * parameters,
	const std::string parameter_name) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | scalar_expr_list_;

    std::size_t num_masks = mask_list_.size();
    p | num_masks;
    if (num_masks != 0) {
      ERROR("Value::pup",
            "Not currently implemented when mask_list_ has 1 or more entries");
    }
  }

  template <class T>
  void evaluate
  (T * values, double t,
   int ndx, int nx, double * x,
   int ndy, int ny, double * y,
   int ndz, int nz, double * z) const throw ();

  double evaluate (double t, double x, double y, double z) const throw ();

  /// returns a string that summarizes contents (for debugging)
  std::string debug_string() const throw();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// List of scalar expressions
  std::vector<ScalarExpr> scalar_expr_list_;

  /// List of masks for each scalar expression
  std::vector< std::shared_ptr<Mask> > mask_list_;

};

#endif /* PROBLEM_VALUE_HPP */

