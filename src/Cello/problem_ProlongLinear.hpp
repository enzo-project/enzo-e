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
  ProlongLinear(CkMigrateMessage *m) : Prolong(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; Prolong::pup(p); }

  /// Prolong coarse Field values to fine values
  virtual void apply
  ( precision_type precision,
    void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
    bool accumulate = false);

  /// Return the name identifying the prolongation operator
  virtual std::string name () const { return "linear"; }

  virtual bool array_sizes_valid (int n3_f[3], int n3_c[3], int * o3) const
  {
    // fine array must have even axes
    const bool is_even =
      ( n3_f[0]%1==0) &&
      ((n3_f[1]%1==0) || n3_f[1]==1) &&
      ((n3_f[2]%1==0) || n3_f[2]==1);

    // fine size must be twice coarse size
    const bool match =
      (n3_f[0] == 2*n3_c[0]) &&
      ((n3_f[1] == 2*n3_c[1]) || (n3_f[1]==1 && n3_c[1]==1)) &&
      ((n3_f[2] == 2*n3_c[2]) || (n3_f[2]==1 && n3_c[2]==1));
      
    return (is_even) && (match);
  }
  
private: // functions

  template <class T>  
  void apply_
  ( T *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const T * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
    bool accumulate = false);

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_PROLONG_LINEAR_HPP */

