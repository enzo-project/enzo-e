// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Prolong.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Problem] Declaration of the Prolong class
///

#ifndef PROBLEM_PROLONG_HPP
#define PROBLEM_PROLONG_HPP

class Prolong 
#ifdef CONFIG_USE_CHARM
  : public PUP::able 
#endif

{

  /// @class    Prolong
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  Prolong() throw();

  /// Destructor
  // ~Prolong() throw();

  // /// Copy constructor
  // Prolong(const Prolong & prolong) throw();

  // /// Assignment operator
  // Prolong & operator= (const Prolong & prolong) throw();

#ifdef CONFIG_USE_CHARM

  /// CHARM++ PUP::able declaration
  PUPable_abstract(Prolong);

  // /// CHARM++ migration constructor for PUP::able
  // Prolong (CkMigrateMessage *m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; PUP::able::pup(p); }

#endif

  /// Prolong comm_block_Ht values to the child block given by (icx,icy,icz)
  virtual void apply 
  (CommBlock        * comm_block_f, 
   const CommBlock  * comm_block_c, 
   const FieldDescr * field_descr,
   int icx, int icy, int icz) = 0;

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_PROLONG_HPP */

