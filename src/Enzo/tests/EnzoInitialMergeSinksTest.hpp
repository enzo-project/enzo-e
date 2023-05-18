// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialMergeSinksTest.hpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2022-01-06
/// @brief    [\ref Enzo] Implementation of test problem for the "merge sinks"
///           method. Creates sink particles with positions, velocities, and
///           masses given in the parameter file 

#ifndef ENZO_ENZO_INITIAL_MERGE_SINKS_TEST_HPP
#define ENZO_ENZO_INITIAL_MERGE_SINKS_TEST_HPP

class EnzoInitialMergeSinksTest : public Initial {

  /// @class    EnzoInitialMergeSinksTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] initial conditions for test problem for "merge sinks"
  ///           method

public: // interface

  /// CHARM++ constructor
  EnzoInitialMergeSinksTest(const EnzoConfig * enzo_config) throw();
 
  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialMergeSinksTest);

  /// CHARM++ migration constructor
  EnzoInitialMergeSinksTest(CkMigrateMessage *m)
    : Initial (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

  // Destructor
  virtual ~EnzoInitialMergeSinksTest() throw() {};

private: // attributes

  // File containing particle data
  std::string particle_data_filename_;

  // Used to store particle data read in from file
  std::vector<enzo_float> mass_data_; 
  std::vector<enzo_float> x_data_; 
  std::vector<enzo_float> y_data_; 
  std::vector<enzo_float> z_data_; 
  std::vector<enzo_float> vx_data_; 
  std::vector<enzo_float> vy_data_; 
  std::vector<enzo_float> vz_data_; 

  // Number of particles
  int n_particles_;

};

#endif /* ENZO_ENZO_INITIAL_MERGE_SINKS_TEST_HPP */

