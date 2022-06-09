// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialGrackleTest.hpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     Tue May 7 2019
/// @brief    [\ref Enzo] Initialization routine for Grackle test problem

#ifndef ENZO_ENZO_INITIAL_ONE_ZONE_FREEFALL_TEST_HPP
#define ENZO_ENZO_INITIAL_ONE_ZONE_FREEFALL_TEST_HPP

class EnzoInitialOneZoneFreefallTest : public Initial {

  /// @class    EnzoInitialGrackleTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for Grackle test problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialOneZoneFreefallTest(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialOneZoneFreefallTest);

  /// CHARM++ migration constructor
  EnzoInitialOneZoneFreefallTest(CkMigrateMessage *m)
     : Initial (m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  (   Block * block, const Hierarchy * hierarchy ) throw();

  // Destructor
  virtual ~EnzoInitialOneZoneFreefallTest() throw() {};

private:

};

#endif /* ENZO_ENZO_INITIAL_ONE_ZONE_FREEFALL_TEST_HPP */
