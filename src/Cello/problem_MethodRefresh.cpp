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
(std::vector< std::string > field_list,
 std::vector< std::string > particle_list,
 int ghost_depth,
 int min_face_rank,
 bool all_fields,
 bool all_particles)
  : Method(),
    field_list_(),
    particle_list_(),
    ghost_depth_(ghost_depth),
    min_face_rank_(min_face_rank),
    all_fields_(all_fields),
    all_particles_(all_particles)
{

  if (field_list.size() > 0) {
    field_list_.resize(field_list.size());
    for (size_t i=0; i<field_list.size(); i++) {
      const int index_field = cello::field_descr()->field_id(field_list[i]);
        field_list_[i] = index_field;
    }
  }

  if (particle_list.size() > 0) {
    particle_list_.resize(particle_list.size());
    for (size_t i=0; i<particle_list.size(); i++) {
      const int index_particle =
        cello::particle_descr()->type_index(particle_list[i]);
      particle_list_[i] = index_particle;
    }
  }

  ASSERT("MethodRefresh()",
         "MethodRefresh requires either field_list or particle_list "
         "to be non-empty",
         (all_fields || all_particles ||
          field_list.size() > 0 || particle_list.size() > 0));

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

