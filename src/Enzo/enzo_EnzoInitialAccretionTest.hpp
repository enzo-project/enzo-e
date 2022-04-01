// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialAccretionTest.hpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2021-03-21
/// @brief    Implementation of EnzoInitialAccretionTest, an initializer for an
///           accretion test problem, which puts an accreting star particle with
///           a given initial position and velocity in a static medium of gas with
///           constant density and internal energy.


#ifndef ENZO_ENZO_INITIAL_ACCRETION_TEST_HPP
#define ENZO_ENZO_INITIAL_ACCRETION_TEST_HPP

class EnzoInitialAccretionTest : public Initial {

  /// @class    EnzoInitialAccretionTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] initial conditions for accretion test problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialAccretionTest
  (int cycle, double time,
   const double star_position[3],
   const double star_velocity[3],
   double star_mass,
   double gas_density,
   double gas_pressure
   ) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialAccretionTest);

  /// CHARM++ migration constructor
  EnzoInitialAccretionTest(CkMigrateMessage *m)
    : Initial (m),
      star_mass_(0.0),
      gas_density_(0.0),
      gas_pressure_(0.0)
  {
    star_position_[0] = 0.0;
    star_position_[1] = 0.0;
    star_position_[2] = 0.0;

    star_velocity_[0] = 0.0;
    star_velocity_[1] = 0.0;
    star_velocity_[2] = 0.0;
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

  // Destructor
  virtual ~EnzoInitialAccretionTest() throw() {};

  private: // attributes
  
  /// Initial position of the star particle
  double star_position_[3];

  /// Initial velocity of the star particle
  double star_velocity_[3];

  /// Initial mass of the star particle
  double star_mass_;

  /// Initial constant density of the gas
  double gas_density_;

  /// Initial constant pressure of the gas
  double gas_pressure_;

};

#endif /* ENZO_ENZO_INITIAL_ACCRETION_TEST_HPP */
