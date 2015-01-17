// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ProlongLinear.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Problem] Declaration of the ProlongLinear class

#ifndef PROBLEM_PROLONG_LINEAR_HPP
#define PROBLEM_PROLONG_LINEAR_HPP

class ProlongLinear : public Prolong 

{

  /// @class    ProlongLinear
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  ProlongLinear() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(ProlongLinear);

  /// CHARM++ migration constructor
  ProlongLinear(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; Prolong::pup(p); }

  /// Prolong fine Field values in the child block (icx,icy,icz) to parent

  virtual int apply
  ( precision_type precision,
    void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3]);

private: // functions

  template <class T>  
  int apply_
  ( T *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const T * values_c, int nd3_c[3], int im3_c[3], int n3_c[3]);

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_PROLONG_LINEAR_HPP */

