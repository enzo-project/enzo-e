// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialMergeSinksTest.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2022-01-06
/// @brief    Implementation of test problem for the "merge sinks" method.
///           Creates sink particles with positions, velocities, and masses
///           given in the file specified by "Initial:merge_sinks_test:particle_data_filename"
///

#include "cello.hpp"
#include "enzo.hpp"
#include <fstream>

//----------------------------------------------------------------------

// #define DEBUG_PERFORMANCE

//----------------------------------------------------------------------

 EnzoInitialMergeSinksTest::EnzoInitialMergeSinksTest
 (const EnzoConfig * enzo_config) throw()
    : Initial (enzo_config->initial_cycle, enzo_config->initial_time)
  {
    // Check sink particle attributes
    cello::particle_descr()->check_particle_attribute("sink","mass");
    cello::particle_descr()->check_particle_attribute("sink","x");
    cello::particle_descr()->check_particle_attribute("sink","y");
    cello::particle_descr()->check_particle_attribute("sink","z");
    cello::particle_descr()->check_particle_attribute("sink","vx");
    cello::particle_descr()->check_particle_attribute("sink","vy");
    cello::particle_descr()->check_particle_attribute("sink","vz");
    cello::particle_descr()->check_particle_attribute("sink","is_copy");
    cello::particle_descr()->check_particle_attribute("sink","id");

    particle_data_filename_ =
    enzo_config->initial_merge_sinks_test_particle_data_filename;

    std::string line;
    std::ifstream inFile(particle_data_filename_);

    n_particles_ = 0;

    while (std::getline(inFile,line)){
      ++n_particles_;
      // Assume each line provides data for one particle
      // Each row must have 7 columns, for mass, x, y, z,
      // vx, vy, vz, respectively
      std::istringstream stream(line);
      enzo_float mass, x, y, z, vx, vy, vz;
      stream >> mass >> x >> y >> z >> vx >> vy >> vz;
      mass_data_.push_back(mass);
      x_data_.push_back(x);
      y_data_.push_back(y);
      z_data_.push_back(z);
      vx_data_.push_back(vx);
      vy_data_.push_back(vy);
      vz_data_.push_back(vz);
    }

    ASSERT("EnzoInitialMergeSinksTest",
	   "Error: No particle data found",
	   n_particles_ != 0);

    return;
  }

void EnzoInitialMergeSinksTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | particle_data_filename_;
  p | mass_data_;
  p | x_data_;
  p | y_data_;
  p | z_data_;
  p | vx_data_;
  p | vy_data_;
  p | vz_data_;
  p | n_particles_;

}

//----------------------------------------------------------------------

void EnzoInitialMergeSinksTest::enforce_block
(Block * block, const Hierarchy * hierarchy ) throw()
{
  // Check if the merge_sinks method is being used
  ASSERT("EnzoInitialMergeSinksTest",
	 "Error: merge_sinks method is required when running with "
	 "the merge_sinks_test initializer.",
	 enzo::problem()->method_exists("merge_sinks"));

  // Check if the pm_update method is being used
  ASSERT("EnzoInitialMergeSinksTest",
	 "Error: pm_update method is required when running with "
	 "the merge_sinks_test initializer.",
	 enzo::problem()->method_exists("pm_update"));
  if (!block->is_leaf()) return;

  ASSERT("EnzoInitialMergeSinksTest",
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

  int it = particle_descr->type_index("sink");

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

  /* Loop through all particles and check if their positions are within
  // bounds of the block. If so, add the particle to this block */

  for (int i = 0 ; i < n_particles_; ++i){

    if (block->check_position_in_block(x_data_[i],
				       y_data_[i],
				       z_data_[i])){

      // Add particle to this block
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

      id[ipp] = i+1; // Want IDs to range from 1 to n_particles_;
      pmass[ipp] = mass_data_[i];
      px[ipp]    = x_data_[i];
      py[ipp]    = y_data_[i];
      pz[ipp]    = z_data_[i];
      pvx[ipp]   = vx_data_[i];
      pvy[ipp]   = vy_data_[i];
      pvz[ipp]   = vz_data_[i];
      is_copy[ipp] = 0;

    } // if particle in block

  } // Loop over particle data
  return;
}
