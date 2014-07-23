// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialTurbulence.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:30:57 UTC 2014
/// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

#ifndef ENZO_ENZO_INITIAL_TURBULENCE_HPP
#define ENZO_ENZO_INITIAL_TURBULENCE_HPP

class EnzoInitialTurbulence : public Initial {

  /// @class    EnzoInitialTurbulence
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialTurbulence() throw() { }
  
  /// Constructor
  EnzoInitialTurbulence(int cycle, double time) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialTurbulence);

  /// CHARM++ migration constructor
  EnzoInitialTurbulence(CkMigrateMessage *m) : Initial (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  (
   CommBlock * block,
   const FieldDescr * field_descr,
   const Hierarchy * hierarchy
   ) throw();

};

#endif /* ENZO_ENZO_INITIAL_TURBULENCE_HPP */

