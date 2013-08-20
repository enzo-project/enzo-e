// See LICENSE_CELLO file for license and copyright information

/// @file     problem_RestrictLinear.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    [\ref Problem] Declaration of the RestrictLinear class

#ifndef PROBLEM_RESTRICT_LINEAR_HPP
#define PROBLEM_RESTRICT_LINEAR_HPP

class RestrictLinear : public Restrict 

{

  /// @class    RestrictLinear
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  RestrictLinear() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(RestrictLinear);

  /// CHARM++ migration constructor
  RestrictLinear(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; PUP::able::pup(p); }

  /// Restrict field_block_ft values to the child block given by (icx,icy,icz)
  int apply 
  ( precision_type precision,
    void *       values_c, int nd3_c[3], int im3_c[3],  int n3_c[3],
    const void * values_f, int nd3_f[3], int im3_f[3],  int n3_f[3]);

private: // functions

  template<class T>
  int apply_
  ( T *       values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
    const T * values_f, int nd3_f[3], int im3_f[3], int n3_f[3]);

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_RESTRICT_LINEAR_HPP */

