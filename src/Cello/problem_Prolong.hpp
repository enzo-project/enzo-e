// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Prolong.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Problem] Declaration of the Prolong class
///

#ifndef PROBLEM_PROLONG_HPP
#define PROBLEM_PROLONG_HPP

class Refresh;
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

  /// CHARM++ migration constructor for PUP::able
  Prolong (CkMigrateMessage *m) :
    PUP::able(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    TRACEPUP;
    PUP::able::pup(p); 
  }

  /// Prolong fine Field values in the child block (icx,icy,icz) to parent

public: // virtual methods
  
  virtual void apply 
  ( precision_type precision,
    void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
    bool accumulate = false) = 0;

  /// Return the name identifying the prolongation operator
  virtual std::string name () const = 0;

  virtual bool array_sizes_valid (int n3_f[3], int n3_c[3], int * o3 = 0) const
  { return true; }

protected: // virtual protected methods

  /// Amount of padding required in coarse region (default 0)
  /// Should only be called by Refresh::coarse_padding()
  friend int Refresh::coarse_padding(const Prolong *) const ;
  virtual int coarse_padding_() const
  { return 0; }

  /// Check whether the size is 

public: // methods
  
protected: // functions


protected: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_PROLONG_HPP */

