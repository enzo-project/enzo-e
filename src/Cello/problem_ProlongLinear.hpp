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
  (      FieldBlock * field_block_f, 
   const FieldBlock * field_block_c, 
   const FieldDescr * field_descr,
   int icx, int icy, int icz);

private: // functions

  template<class T>
  void interpolate_(T       * values_f,
		    const T * values_c,
		    int ndx, int ndy, int ndz,
		    int ixm, int iym, int izm,
		    int nx, int ny, int nz,
		    int gx, int gy, int gz);

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_PROLONG_LINEAR_HPP */

