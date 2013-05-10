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

#ifdef CONFIG_USE_CHARM

  /// CHARM++ PUP::able declaration
  PUPable_decl(RestrictLinear);

  /// CHARM++ migration constructor
  RestrictLinear(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; PUP::able::pup(p); }

#endif

  /// Restrict comm_block_ft values to the child block given by (icx,icy,icz)
  virtual void apply 
  (CommBlock        * comm_block_c, 
   const CommBlock  * comm_block_f, 
   const FieldDescr * field_descr,
   int icx, int icy, int icz);

private: // functions

  template<class T>
  void interpolate_(T       * values_c,
		    const T * values_f,
		    int ndx, int ndy, int ndz,
		    int ixm, int iym, int izm,
		    int nx, int ny, int nz,
		    int gx, int gy, int gz);

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_RESTRICT_LINEAR_HPP */

