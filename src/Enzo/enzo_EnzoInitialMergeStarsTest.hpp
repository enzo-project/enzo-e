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

#endif /* ENZO_ENZO_INITIAL_MERGE_STARS_TEST_HPP */

