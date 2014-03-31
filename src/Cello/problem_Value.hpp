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
  Value() throw();

  /// Destructor
  ~Value() throw();

  /// Copy constructor
  Value(const Value & Value) throw();

  /// Assignment operator
  Value & operator= (const Value & Value) throw();

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

  // NOTE: change pup() function whenever attributes change

  /// List of scalar expressions
  std::vector<ScalarExpr *> scalar_expr_;

  /// List of masks for each scalar expression
  std::vector<Mask *> mask_;

};

#endif /* PROBLEM_VALUE_HPP */

