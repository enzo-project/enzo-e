// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Prolong.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Problem] Declaration of the Prolong class
///

#ifndef PROBLEM_PROLONG_HPP
#define PROBLEM_PROLONG_HPP

class Prolong : public PUP::able 

{

  /// @class    Prolong
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  Prolong() throw();

  /// CHARM++ PUP::able declaration
  PUPable_abstract(Prolong);

  // /// CHARM++ migration constructor for PUP::able
  // Prolong (CkMigrateMessage *m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    TRACEPUP;
    PUP::able::pup(p); 
    p | monotonic_;
    p | positive_;
  }

  /// Prolong fine Field values in the child block (icx,icy,icz) to parent

  virtual int apply 
  ( precision_type precision,
    void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3]) = 0;

  /// Set whether interpolation should be monotonic
  void set_monotonic (bool monotonic) 
  { monotonic_ = monotonic; }

  /// Return monotonicity setting
  bool monotonic (bool monotonic) const
  { return monotonic_; }

  /// Set whether interpolation should be positive
  void set_positive (bool positive) 
  { positive_ = positive; }

  /// Return positivity setting
  bool positive (bool positive) const
  { return positive_; }

protected: // functions


protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Whether interpolation should be monotonic
  bool monotonic_;

  /// Whether interpolation should be positive
  bool positive_;

};

#endif /* PROBLEM_PROLONG_HPP */

