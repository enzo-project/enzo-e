// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jul  8 12:22:25 PDT 2014
/// @brief    [\ref Enzo] Initialization routine for Grackle test problem

#ifndef ENZO_ENZO_INITIAL_GRACKLE_TEST_HPP
#define ENZO_ENZO_INITIAL_GRACKLE_TEST_HPP

class EnzoInitialGrackleTest : public Initial {

  /// @class    EnzoInitialGrackleTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for Grackle test problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialGrackleTest() throw() { }
  
  /// Constructor
  EnzoInitialGrackleTest(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialGrackleTest);

  /// CHARM++ migration constructor
  EnzoInitialGrackleTest(CkMigrateMessage *m) : Initial (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  (
   Block * block,
   const FieldDescr * field_descr,
   const Hierarchy * hierarchy
   ) throw();

private:

#ifdef CONFIG_USE_GRACKLE

  const code_units      * units_;
  const chemistry_data  * chemistry_;

#endif
  
};

#endif /* ENZO_ENZO_INITIAL_GRACKLE_TEST_HPP */

