// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ProlongElliptic.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Problem] Declaration of the ProlongElliptic class

#ifndef PROBLEM_PROLONG_ELLIPTIC_HPP
#define PROBLEM_PROLONG_ELLIPTIC_HPP

class ProlongElliptic : public Prolong 

{

  /// @class    ProlongElliptic
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  ProlongElliptic() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(ProlongElliptic);

  /// CHARM++ migration constructor
  ProlongElliptic(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; Prolong::pup(p); }

  /// Prolong comm_block_Ht values to the child block given by (icx,icy,icz)
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

#endif /* PROBLEM_PROLONG_ELLIPTIC_HPP */

