// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ProlongInject.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-01-10
/// @brief    [\ref Problem] Declaration of the ProlongInject class

#ifndef PROBLEM_PROLONG_INJECT_HPP
#define PROBLEM_PROLONG_INJECT_HPP

class ProlongInject : public Prolong 

{

  /// @class    ProlongInject
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  ProlongInject() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(ProlongInject);

  /// CHARM++ migration constructor
  ProlongInject(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; Prolong::pup(p); }

  /// Prolong fine Field values in the child block (icx,icy,icz) to parent

  virtual int apply
  ( precision_type precision,
    void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3]);

  /// Return the name identifying the prolongation operator
  virtual std::string name () const { return "inject"; }

private: // functions

  template <class T>  
  int apply_
  ( T *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const T * values_c, int nd3_c[3], int im3_c[3], int n3_c[3]);

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_PROLONG_INJECT_HPP */

