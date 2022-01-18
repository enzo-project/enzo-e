// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialMergeStarsTest.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2022-01-06
/// @brief    Implementation of test problem for the "merge stars" method.
///           Sets up two star particles with positions, velocities, and masses
///           given in the parameter file
///

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

// #define DEBUG_PERFORMANCE

//----------------------------------------------------------------------

 EnzoInitialMergeStarsTest::EnzoInitialMergeStarsTest 
 (const EnzoConfig * enzo_config) throw()
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
  
  EnzoSimulation * enzo_simulation = enzo::simulation();
  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);

  // Initialize particle data

  ParticleDescr * particle_descr = cello::particle_descr();
  Particle particle = block->data()->particle();

  int it = particle_descr->type_index("star");

  int ia_m = particle.attribute_index (it, "mass");
  int ia_x = particle.attribute_index (it, "x");
  int ia_y = particle.attribute_index (it, "y");
  int ia_z = particle.attribute_index (it, "z");
  int ia_vx = particle.attribute_index (it, "vx");
  int ia_vy = particle.attribute_index (it, "vy");
  int ia_vz = particle.attribute_index (it, "vz");
  int ia_copy  = particle.attribute_index (it, "is_copy");
  int ia_id   = particle.attribute_index (it, "id");

  int ib  = 0; // batch counter
  int ipp = 0; // particle counter

  // this will point to the particular value in the
  // particle attribute array
  enzo_float * pmass = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
  int64_t * is_copy = 0;
  int64_t * id = 0;

  // for both particles, check if their positions are within the 
  // bounds of the block. If so, add the particle to this block

  if (block->check_position_in_block(pos_1_[0],pos_1_[1],pos_1_[2])){
  
    // Add particle 1 to this block
    int new_particle_index = particle.insert_particles(it, 1);
    // Add particle to simulation
    enzo_simulation->data_insert_particles(1);
    
    particle.index(new_particle_index,&ib,&ipp);

    id   = (int64_t *) particle.attribute_array(it, ia_id, ib);
    pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
    px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
    pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
    pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
    pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
    is_copy   = (int64_t *) particle.attribute_array(it, ia_copy, ib);

    id[ipp] = 0;
    pmass[ipp] = mass_1_;
    px[ipp]    = pos_1_[0];
    py[ipp]    = pos_1_[1];
    pz[ipp]    = pos_1_[2];
    pvx[ipp]   = vel_1_[0];
    pvy[ipp]   = vel_1_[1];
    pvz[ipp]   = vel_1_[2];
    is_copy[ipp] = 1;
  
  
  } // if particle 1 in block

  if (block->check_position_in_block(pos_2_[0],pos_2_[1],pos_2_[2])){
  
    // Add particle 1 to this block
    int new_particle_index = particle.insert_particles(it, 1);
    // Add particle to simulation
    enzo_simulation->data_insert_particles(1);
    
    particle.index(new_particle_index,&ib,&ipp);

    id   = (int64_t *) particle.attribute_array(it, ia_id, ib);
    pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
    px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
    pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
    pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
    pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
    is_copy   = (int64_t *) particle.attribute_array(it, ia_copy, ib);

    id[ipp] = 0;
    pmass[ipp] = mass_2_;
    px[ipp]    = pos_2_[0];
    py[ipp]    = pos_2_[1];
    pz[ipp]    = pos_2_[2];
    pvx[ipp]   = vel_2_[0];
    pvy[ipp]   = vel_2_[1];
    pvz[ipp]   = vel_2_[2];
    is_copy[ipp] = 1;
  
  
  } // if particle 2 in block

  return;
}


