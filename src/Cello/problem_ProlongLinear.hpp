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

#ifdef CONFIG_USE_CHARM

  /// CHARM++ PUP::able declaration
  PUPable_decl(ProlongLinear);

  /// CHARM++ migration constructor
  ProlongLinear(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; PUP::able::pup(p); }

#endif

  /// Prolong comm_block_Ht values to the child block given by (icx,icy,icz)
  virtual void apply 
  ( precision_type precision,
     void * values_f,
    int ndx_f, int ndy_f, int ndf_z, 
    int nx_f,  int ny_f,  int nz_f,
    const void * values_c,
    int ndx_c, int ndy_c, int ndc_z, 
    int nx_c,  int ny_c,  int nz_c);

private: // functions

  template <class T>  
  void apply_
  ( T * values_f,
    int ndx_f, int ndy_f, int ndf_z, 
    int nx_f,  int ny_f,  int nz_f,
    const T * values_c,
    int ndx_c, int ndy_c, int ndc_z, 
    int nx_c,  int ny_c,  int nz_c);

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_PROLONG_LINEAR_HPP */

