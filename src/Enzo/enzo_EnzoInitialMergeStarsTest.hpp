// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialMergeStarsTest.hpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2022-01-06
/// @brief    [\ref Enzo] Implementation of test problem for the "merge stars"
///           method. Sets up two star particles with positions, velocities, and
///           masses given in the parameter file 

#ifndef ENZO_ENZO_INITIAL_MERGE_STARS_TEST_HPP
#define ENZO_ENZO_INITIAL_MERGE_STARS_TEST_HPP

class EnzoInitialMergeStarsTest : public Initial {

  /// @class    EnzoInitialMergeStarsTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] initial conditions for test problem for "merge stars"
  ///           method

public: // interface

  /// CHARM++ constructor
  EnzoInitialMergeStarsTest(const EnzoConfig * enzo_config) throw();
 
  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialMergeStarsTest);

  /// CHARM++ migration constructor
  EnzoInitialMergeStarsTest(CkMigrateMessage *m)
    : Initial (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

  // Destructor
  virtual ~EnzoInitialMergeStarsTest() throw() {};

private: // attributes

  // Position of first particle
  double pos_1_[3];

  // Position of second particle
  double pos_2_[3];

  // Velocity of first particle
  double vel_1_[3];

  // Velocity of second particle
  double vel_2_[3];

  // Mass of first particle
  double mass_1_;

  // Mass of second particle
  double mass_2_;

};

#endif /* ENZO_ENZO_INITIAL_MERGE_STARS_TEST_HPP */

