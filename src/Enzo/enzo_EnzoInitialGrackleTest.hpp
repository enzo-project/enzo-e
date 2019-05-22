// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.hpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     Tue May 7 2019
/// @brief    [\ref Enzo] Initialization routine for Grackle test problem

#ifndef ENZO_ENZO_INITIAL_GRACKLE_TEST_HPP
#define ENZO_ENZO_INITIAL_GRACKLE_TEST_HPP

class EnzoInitialGrackleTest : public Initial {

  /// @class    EnzoInitialGrackleTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for Grackle test problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialGrackleTest(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialGrackleTest);

  /// CHARM++ migration constructor
  EnzoInitialGrackleTest(CkMigrateMessage *m)
     : Initial (m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  (   Block * block, const Hierarchy * hierarchy ) throw();

  // Destructor
  virtual ~EnzoInitialGrackleTest() throw() {};

private:

};

#endif /* ENZO_ENZO_INITIAL_GRACKLE_TEST_HPP */
