// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialImplosion2.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

#ifndef ENZO_ENZO_INITIAL_IMPLOSION2_HPP
#define ENZO_ENZO_INITIAL_IMPLOSION2_HPP

class EnzoInitialImplosion2 : public Initial {

  /// @class    EnzoInitialImplosion2
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialImplosion2() throw() { }
  
  /// Constructor
  EnzoInitialImplosion2(int cycle, double time) throw();

#ifdef CONFIG_USE_CHARM

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialImplosion2);

  /// CHARM++ migration constructor
  EnzoInitialImplosion2(CkMigrateMessage *m) : Initial (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif

  /// Initialize the block

  virtual void enforce 
  (
   Block * block,
   const FieldDescr * field_descr,
   const Hierarchy * hierarchy
   ) throw();

};

#endif /* ENZO_ENZO_INITIAL_IMPLOSION2_HPP */

