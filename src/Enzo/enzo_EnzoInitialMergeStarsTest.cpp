// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialMergeStarsTest.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2022-01-06
/// @brief    Implementation of test problem for the "merge stars" method.
///           Sets up two scodestar particles with positions, velocities, and masses
///           given in the parameter file
///

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

// #define DEBUG_PERFORMANCE

//----------------------------------------------------------------------

 EnzoInitialMergeStarsTest::EnzoInitialMergeStarsTest (const EnzoConfig * enzo_config) throw()
    : Initial (enzo_config->initial_cycle, enzo_config->initial_time)
  {
    pos_1_[0] = enzo_config->initial_merge_stars_test_pos_1[0];
    pos_1_[1] = enzo_config->initial_merge_stars_test_pos_1[1];
    pos_1_[2] = enzo_config->initial_merge_stars_test_pos_1[2];
    
    pos_2_[0] = enzo_config->initial_merge_stars_test_pos_2[0];
    pos_2_[1] = enzo_config->initial_merge_stars_test_pos_2[1];
    pos_2_[2] = enzo_config->initial_merge_stars_test_pos_2[2];
    
    vel_1_[0] = enzo_config->initial_merge_stars_test_vel_1[0];
    vel_1_[1] = enzo_config->initial_merge_stars_test_vel_1[1];
    vel_1_[2] = enzo_config->initial_merge_stars_test_vel_1[2];
    
    vel_2_[0] = enzo_config->initial_merge_stars_test_vel_2[0];
    vel_2_[1] = enzo_config->initial_merge_stars_test_vel_2[1];
    vel_2_[2] = enzo_config->initial_merge_stars_test_vel_2[2];
  
    mass_1_ = enzo_config->initial_merge_stars_test_mass_1;
    mass_2_ = enzo_config->initial_merge_stars_test_mass_2;
  }

void EnzoInitialMergeStarsTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,pos_1_,3);
  PUParray(p,pos_2_,3);
  PUParray(p,vel_1_,3);
  PUParray(p,vel_2_,3);
  p | mass_1_;
  p | mass_2_;
  
}

//----------------------------------------------------------------------

void EnzoInitialMergeStarsTest::enforce_block
(Block * block, const Hierarchy * hierarchy ) throw()
{

  if (!block->is_leaf()) return;

  ASSERT("EnzoInitialMergeStarsTest",
  	 "Block does not exist",
  	 block != NULL);
  
  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);


  // Initialize particles

  Particle particle = block->data()->particle();
  
}

