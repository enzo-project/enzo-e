// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodRefresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-03-09
/// @brief    Implementation of the refresh "method"

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

//----------------------------------------------------------------------

MethodRefresh::MethodRefresh
(std::vector<int> field_list,
 std::vector<int> particle_list,
 int ghost_depth,
 int min_face_rank,
 bool all_fields,
 bool all_particles)
  : Method(),
    field_list_(field_list),
    particle_list_(particle_list),
    ghost_depth_(ghost_depth),
    min_face_rank_(min_face_rank),
    all_particles_(all_particles),
    all_fields_(all_fields)
{

  Refresh * refresh = cello::refresh(ir_post_);

  if (ghost_depth > 0) refresh->set_ghost_depth (ghost_depth);
  refresh->set_min_face_rank(min_face_rank);

  // add fields to refresh
  if (all_fields_) {
    refresh->add_all_fields();
  } else {
    const int nf=field_list_.size();
    for (int i_f=0; i_f<nf; i_f++) {
      cello::refresh(ir_post_)->add_field(field_list_[i_f]);
    }
  }

  // add particles to refresh
  if (all_particles_) {
    refresh->add_all_particles();
  } else {
    const int nf=particle_list_.size();
    for (int i_f=0; i_f<nf; i_f++) {
      cello::refresh(ir_post_)->add_particle(particle_list_[i_f]);
    }
  }
}

//----------------------------------------------------------------------

void MethodRefresh::compute ( Block * block) throw()
{
  block->compute_done(); 
}

